"""Exploratory ECM NEON timing screen for a frozen Stage A cell."""

import argparse
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path

from ....ecm.neon_probe import ARMS, ENTRIES, WIDTHS, NEONProbe, file_sha

IDENTITY = "b4a-ecm-neon-feasibility-v3"
CELL_ID = "ecm-neon-feasibility-v3"
COLUMNS = (
    "phase",
    "width",
    "operation",
    "arm",
    "cycle",
    "repeats",
    "offset",
    "wall",
    "cpu",
    "pairs",
    "checksum",
    "validated",
    "status",
    "qos_set",
    "qos_get",
    "requested_qos",
    "observed_qos",
    "relative_priority",
    "thread_id",
    "thread_id_status",
    "independent_operations",
    "since_entry",
)


def screen_cell():
    return dict(
        id=CELL_ID,
        arm="ECM-NEON-screen",
        endpoint="neon-screen",
        mode="exploratory",
        wall=180,
        cpu=360,
        shutdown=2,
        rss=4 * 2**30,
        output=8 * 2**20,
        processes=2,
        threads=8,
        seed=2026110407,
        warmup=3,
        samples=9,
        minimum_batch=0.01,
        maximum_batch=0.1,
        min_repeats=64,
        max_repeats=65536,
        max_rows=16000,
        raw_bytes=6 * 2**20,
    )


def check_screen_cell(cell):
    if cell != screen_cell():
        raise ValueError(
            "changed NEON screen needs a new prospective protocol"
        )


def panel(width, seed):
    """Fixed full-width modulus and varied canonical residues, not a corpus."""
    n = 2**width - 159
    generator = random.Random(seed + width)
    # Boundary/carry fixtures belong to correctness admission. Timed operands
    # are uniformly drawn canonical residues, without forcing high bits.
    rows = [
        tuple(generator.randrange(n) for _ in range(4)) for _ in range(ENTRIES)
    ]
    return n, rows


def summarize(rows, seed):
    """Separate matched arithmetic, hybrid and layout costs; no adoption."""
    generator, output = random.Random(seed), []
    for width in WIDTHS:
        for operation in ("multiply", "square"):
            selected = [
                r
                for r in rows
                if r["width"] == width and r["operation"] == operation
            ]
            by_arm = {
                arm: sorted(
                    (r for r in selected if r["arm"] == arm),
                    key=lambda r: r["cycle"],
                )
                for arm in ARMS
            }
            summaries = {}
            for arm, observations in by_arm.items():
                if len(observations) != 9 or {
                    r["cycle"] for r in observations
                } != set(range(9)):
                    raise ValueError(
                        "nine complete matched batch observations required"
                    )
                if any(
                    r["independent_operations"]
                    != (0 if arm == ARMS[6] else 2 * r["pairs"])
                    for r in observations
                ):
                    raise ValueError("independent arithmetic count changed")
                values = [r["wall"] / (2 * r["pairs"]) for r in observations]
                median = statistics.median(values)
                summaries[arm] = dict(
                    median_seconds_per_result=median,
                    min=min(values),
                    max=max(values),
                    mad=statistics.median(abs(x - median) for x in values),
                )
            for group in (ARMS[:3], ARMS[3:]):
                for cycle in range(9):
                    if (
                        len({by_arm[arm][cycle]["pairs"] for arm in group})
                        != 1
                        or len({by_arm[arm][cycle]["offset"] for arm in group})
                        != 1
                    ):
                        raise ValueError(
                            "unmatched operation counts or input order"
                        )
            ratios = [
                new["wall"] / old["wall"]
                for new, old in zip(by_arm[ARMS[5]], by_arm[ARMS[3]])
            ]
            boots = sorted(
                statistics.median(generator.choices(ratios, k=9))
                for _ in range(2000)
            )
            qos = {
                (
                    r["qos_set"],
                    r["qos_get"],
                    r["requested_qos"],
                    r["observed_qos"],
                    r["relative_priority"],
                    r["thread_id"],
                    r["thread_id_status"],
                )
                for r in selected
            }
            qos_compatible = len(qos) == 1 and all(
                r["qos_set"] == r["qos_get"] == r["thread_id_status"] == 0
                and r["requested_qos"] == r["observed_qos"] == 0x21
                and r["relative_priority"] == 0
                for r in selected
            )

            def ratio(a, b):
                return (
                    summaries[a]["median_seconds_per_result"]
                    / summaries[b]["median_seconds_per_result"]
                )

            output.append(
                dict(
                    width=width,
                    operation=operation,
                    arms=summaries,
                    arithmetic_neon_over_mpn=ratio(ARMS[2], ARMS[0]),
                    arithmetic_neon_over_same_layout_scalar=ratio(
                        ARMS[2], ARMS[1]
                    ),
                    complete_neon_over_production=ratio(ARMS[5], ARMS[3]),
                    complete_neon_over_same_layout_scalar=ratio(
                        ARMS[5], ARMS[4]
                    ),
                    complete_scalar_layout_and_algorithm_over_production=ratio(
                        ARMS[4], ARMS[3]
                    ),
                    layout_only_diagnostic_over_production=ratio(
                        ARMS[6], ARMS[3]
                    ),
                    paired_complete_ratio_interval=[boots[49], boots[1949]],
                    qos_compatible=qos_compatible,
                    supports_further_integration=qos_compatible
                    and boots[1949] < 1,
                )
            )
    return output


def validate_report(data, cell):
    check_screen_cell(cell)
    if (
        data.get("identity") != IDENTITY
        or type(data.get("complete")) is not bool
    ):
        raise ValueError("invalid NEON screen output")
    if not data["complete"]:
        return
    rows = data["measurements"]
    if len(rows) != len(WIDTHS) * 2 * len(ARMS) * 9:
        raise ValueError("incomplete NEON measurement grid")
    for row in rows:
        for name in ("wall", "cpu"):
            if (
                type(row[name]) not in (int, float)
                or not math.isfinite(row[name])
                or row[name] <= 0
            ):
                raise ValueError("invalid native timing")
        if (
            row["status"] != 0
            or row["validated"] is not True
            or row["independent_operations"]
            != (0 if row["arm"] == ARMS[6] else row["repeats"] * ENTRIES * 2)
        ):
            raise ValueError("unvalidated or incomplete native batch")
    warm = data["warmup"]
    expected_warm = {
        (width, operation, arm)
        for width in WIDTHS
        for operation in ("multiply", "square")
        for arm in ARMS
    }
    if (
        len(warm) != len(expected_warm)
        or {(r["width"], r["operation"], r["arm"]) for r in warm}
        != expected_warm
        or any(
            type(r["validated_driver_seconds"]) not in (int, float)
            or not math.isfinite(r["validated_driver_seconds"])
            or r["validated_driver_seconds"] < 3
            for r in warm
        )
    ):
        raise ValueError("required validated warmup missing")
    if data["summary"] != summarize(rows, cell["seed"]):
        raise ValueError("NEON summary does not reconstruct from raw batches")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy Python 3.11 required")
    request = json.loads(args.request.read_text())
    check_screen_cell(request["cell"])
    if (
        request.get("stage") != "A"
        or request.get("admission_passed") is not True
    ):
        raise ValueError("frozen Stage A correctness admission required")
    cell = request["cell"]
    started = request["entry"]
    raw_path = Path(request["raw_path"])
    measured, warmup, raw_rows, raw_bytes = [], [], 0, 0
    result = dict(
        identity=IDENTITY,
        complete=False,
        measurements=measured,
        warmup=warmup,
        native_threads=1,
        placement="unobserved; QoS is not P/E pinning",
    )

    def poll():
        if time.monotonic() - started >= cell["wall"] - cell["shutdown"] - 3:
            raise RuntimeError("screen_lifetime_limit")

    with raw_path.open("xb") as raw:
        header = dict(
            columns=COLUMNS,
            arms=ARMS,
            identity=IDENTITY,
            operation_names=["multiply", "square"],
            qos_request=(
                "pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0)"
            ),
            qos_readback="requested class; overrides and placement unobserved",
        )
        raw.write((json.dumps(header, separators=(",", ":")) + "\n").encode())

        def observe(probe, arm, operation, repeats, offset, phase, cycle):
            nonlocal raw_rows, raw_bytes
            poll()
            before = time.monotonic()
            row = probe.batch(arm, operation, repeats, offset)
            row.update(
                width=probe.width,
                operation=operation,
                arm=arm,
                phase=phase,
                cycle=cycle,
                since_entry=before - started,
            )
            packed = dict(
                row, operation=(operation == "square"), arm=ARMS.index(arm)
            )
            encoded = (
                json.dumps(
                    [packed[key] for key in COLUMNS],
                    separators=(",", ":"),
                    allow_nan=False,
                )
                + "\n"
            ).encode()
            raw_rows += 1
            raw_bytes += len(encoded)
            if raw_rows > cell["max_rows"] or raw_bytes > cell["raw_bytes"]:
                raise RuntimeError("recorder_exhausted")
            raw.write(encoded)
            raw.flush()
            if not row["validated"] or row["status"]:
                raise RuntimeError("invalid_or_truncated_native_batch")
            poll()
            return row, time.monotonic() - before

        try:
            cases = [
                (width, operation)
                for width in WIDTHS
                for operation in ("multiply", "square")
            ]
            random.Random(cell["seed"]).shuffle(cases)
            for width, operation in cases:
                n, values = panel(width, cell["seed"])
                probe = NEONProbe(width, n)
                try:
                    for index, inputs in enumerate(values):
                        probe.set(index, inputs)
                    # Untimed admission on this exact panel before any clocked
                    # calibration. Each arm is also checked internally vs GMP.
                    for arm in ARMS:
                        for index in range(ENTRIES):
                            probe.once(arm, operation, index)
                    chosen = {}
                    generator = random.Random(
                        cell["seed"] + width + (operation == "square")
                    )
                    for group in (ARMS[:3], ARMS[3:]):
                        repeats = cell["min_repeats"]
                        for cycle in range(11):
                            order = list(group)
                            generator.shuffle(order)
                            rows = [
                                observe(
                                    probe,
                                    arm,
                                    operation,
                                    repeats,
                                    cycle % ENTRIES,
                                    "calibrate",
                                    cycle,
                                )[0]
                                for arm in order
                            ]
                            fastest, slowest = (
                                min(r["wall"] for r in rows),
                                max(r["wall"] for r in rows),
                            )
                            if slowest > cell["maximum_batch"]:
                                raise RuntimeError(
                                    "matched_batch_spread_infeasible"
                                )
                            if fastest >= cell["minimum_batch"]:
                                chosen.update({arm: repeats for arm in group})
                                break
                            repeats *= 2
                            if repeats > cell["max_repeats"]:
                                raise RuntimeError(
                                    "calibration_repetition_limit"
                                )
                        else:
                            raise RuntimeError("calibration_call_limit")
                    accumulated = {arm: 0.0 for arm in ARMS}
                    calls = {arm: 0 for arm in ARMS}
                    while min(accumulated.values()) < cell["warmup"]:
                        order = [
                            arm
                            for arm in ARMS
                            if accumulated[arm] < cell["warmup"]
                        ]
                        generator.shuffle(order)
                        for arm in order:
                            calls[arm] += 1
                            if calls[arm] > 512:
                                raise RuntimeError("warmup_call_limit")
                            _, duration = observe(
                                probe,
                                arm,
                                operation,
                                chosen[arm],
                                calls[arm] % ENTRIES,
                                "warmup",
                                calls[arm],
                            )
                            accumulated[arm] += duration
                    warmup.extend(
                        dict(
                            width=width,
                            operation=operation,
                            arm=arm,
                            validated_driver_seconds=accumulated[arm],
                            calls=calls[arm],
                        )
                        for arm in ARMS
                    )
                    for group in (ARMS[:3], ARMS[3:]):
                        repeats = chosen[group[0]]
                        for cycle in range(11):
                            order = list(group)
                            generator.shuffle(order)
                            rows = [
                                observe(
                                    probe,
                                    arm,
                                    operation,
                                    repeats,
                                    cycle % ENTRIES,
                                    "calibrate-warm",
                                    cycle,
                                )[0]
                                for arm in order
                            ]
                            if (
                                max(r["wall"] for r in rows)
                                > cell["maximum_batch"]
                            ):
                                raise RuntimeError(
                                    "warm_batch_spread_infeasible"
                                )
                            if (
                                min(r["wall"] for r in rows)
                                >= cell["minimum_batch"]
                            ):
                                chosen.update({arm: repeats for arm in group})
                                break
                            repeats *= 2
                            if repeats > cell["max_repeats"]:
                                raise RuntimeError("warm_calibration_limit")
                        else:
                            raise RuntimeError("warm_calibration_call_limit")
                    base_order = list(ARMS)
                    generator.shuffle(base_order)
                    for cycle in range(cell["samples"]):
                        order = (
                            base_order[cycle % len(ARMS) :]
                            + base_order[: cycle % len(ARMS)]
                        )
                        if cycle % 2:
                            order.reverse()
                        for arm in order:
                            row, _ = observe(
                                probe,
                                arm,
                                operation,
                                chosen[arm],
                                (17 * cycle + width) % ENTRIES,
                                "measure",
                                cycle,
                            )
                            measured.append(row)
                finally:
                    probe.close()
            result.update(
                complete=True, summary=summarize(measured, cell["seed"])
            )
        except BaseException as error:
            result.update(
                complete=False,
                failure=type(error).__name__ + ": " + str(error),
            )
        finally:
            raw.flush()
    result.update(
        raw_sha256=file_sha(raw_path),
        raw_rows=raw_rows,
        elapsed=time.monotonic() - started,
    )
    print(json.dumps(result, separators=(",", ":"), allow_nan=False))


if __name__ == "__main__":
    main()
