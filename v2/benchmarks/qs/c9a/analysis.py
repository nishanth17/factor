"""Modulus-level uncertainty with explicit censoring and drift failures."""

import math
import random
import statistics
import time
from collections import defaultdict


def interval(values):
    ordered = sorted(values)
    return [
        ordered[int(0.025 * (len(ordered) - 1))],
        ordered[int(0.975 * (len(ordered) - 1))],
    ]


def summarize(results, *, nominee="C1", samples=4096):
    """Bootstrap independent moduli, averaging repeated seeds within each."""
    if type(samples) is not int or not 100 <= samples <= 10000:
        raise ValueError("bootstrap allowance is outside its finite plan")
    started, cpu = time.monotonic(), time.process_time()
    groups = defaultdict(lambda: defaultdict(dict))
    for result in results:
        if result["arm"] not in ("C0", nominee):
            continue
        key = (result["digits"], result.get("uneven", False))
        identity = (result["fixture"], result["seed"], result["arm"])
        if identity in groups[key][result["fixture"]]:
            raise ValueError("duplicate outcome would inflate independence")
        groups[key][result["fixture"]][identity] = result
    output, work = [], len(results)
    rng = random.Random(202611030903)
    for (digits, uneven), moduli in sorted(groups.items()):
        completion, costs, ratios, invalid, missing = [], [], [], 0, 0
        for fixture, cells in moduli.items():
            seeds = sorted({key[1] for key in cells})
            differences, cost, times = [], [], []
            for seed in seeds:
                left = cells.get((fixture, seed, "C0"))
                right = cells.get((fixture, seed, nominee))
                if left is None or right is None:
                    missing += 1
                    continue
                differences.append(
                    int(right["complete"]) - int(left["complete"])
                )
                valid = all(
                    r.get("timing_valid", False)
                    and r.get("observation", "clean") == "clean"
                    for r in (left, right)
                )
                if not valid:
                    invalid += 1
                    continue
                cost.append(right["capped_seconds"] / left["capped_seconds"])
                if left["complete"] and right["complete"]:
                    times.append(math.log(right["seconds"] / left["seconds"]))
            if differences:
                completion.append(statistics.mean(differences))
            if len(cost) == len(seeds):
                costs.append(statistics.mean(cost))
            if len(times) == len(seeds):
                ratios.append(statistics.mean(times))
        record = dict(
            digits=digits,
            uneven=uneven,
            independent_moduli=len(moduli),
            paired_moduli=len(completion),
            invalid_timing_pairs=invalid,
            missing_pairs=missing,
            complete_cost_moduli=len(ratios),
            completion_effect=None
            if not completion
            else statistics.mean(completion),
            observed_completion_regression_over_5pp=bool(
                completion and statistics.mean(completion) < -0.05
            ),
            balanced_60_70_completion_regression=bool(
                digits in (60, 70)
                and not uneven
                and completion
                and statistics.mean(completion) < 0
            ),
        )
        for name, values in (
            ("completion_difference", completion),
            ("capped_cost_ratio", costs),
            ("log_time_ratio", ratios),
        ):
            if not values:
                record[name + "_ci95"] = None
                continue
            boot = []
            for _ in range(samples):
                if (
                    max(time.monotonic() - started, time.process_time() - cpu)
                    > 110
                ):
                    raise RuntimeError("analysis wall/CPU allowance exhausted")
                boot.append(
                    statistics.mean(rng.choice(values) for _ in values)
                )
                work += 16 * len(values)
            bounds = interval(boot)
            record[name + "_ci95"] = bounds
            if name == "log_time_ratio":
                record["complete_time_ratio_ci95"] = [
                    math.exp(x) for x in bounds
                ]
        # A degenerate bootstrap interval is not proof of population
        # noninferiority. Even zero regressions among n independent moduli
        # has this one-sided 95% upper bound on unseen regression frequency.
        n = len(completion)
        record["zero_event_regression_upper95"] = (
            1 - 0.05 ** (1 / n)
            if n and all(v >= 0 for v in completion)
            else None
        )
        record["confirmation_precision"] = (
            "limited" if n < 9 or uneven else "review CI and censoring"
        )
        output.append(record)
    return dict(
        classes=output,
        work=work,
        bootstrap_samples=samples,
        inference="seeds nested in modulus; no automatic adoption",
        invalid_policy="retain all outcomes; count excluded invalid times",
    )


def observer_summary(results):
    """Nine identical-work overhead pairs and bounded stream capacity facts."""
    if len(results) != 18:
        raise ValueError("observer check requires all nine predeclared pairs")
    ratios, cpu_ratios, byte_counts = [], [], []
    for index in range(9):
        pair = results[2 * index : 2 * index + 2]
        clean = next(r for r in pair if r["observation"] == "clean")
        capture = next(r for r in pair if r["observation"] == "capture")

        def signature(result):
            return (
                result["divisor"],
                result["reason"],
                result["next_position"],
                result["measured_work"],
                result["stats"].get("engine", {}).get("collector_counts"),
            )

        if signature(clean) != signature(capture):
            raise ValueError("observer compared different physical work")
        if (
            not all(r["timing_valid"] for r in pair)
            or capture["recorder"]["exhausted"]
        ):
            raise ValueError(
                "observer drift/resource screen failed; no automatic retry"
            )
        ratios.append(capture["seconds"] / clean["seconds"])
        cpu_ratios.append(capture["cpu_seconds"] / clean["cpu_seconds"])
        byte_counts.append(capture["recorder"]["bytes"])
    return dict(
        work=18 * 1024,
        overhead_ratio=statistics.median(ratios),
        cpu_overhead_ratio=statistics.median(cpu_ratios),
        ratios=ratios,
        cpu_ratios=cpu_ratios,
        stream_bytes=byte_counts,
        capacity_projection_required=True,
        scope="one block; mature capture capacity still needs justification",
    )


def intake_summary(results):
    """Expose physical coverage, rejection costs and useful-cycle proxies."""
    cells, coverage = [], defaultdict(set)
    for row in results:
        if (
            row.get("observation") != "clean"
            or row.get("cap_seconds", 0) <= 10
        ):
            continue
        engine = row.get("stats", {}).get("engine", {})
        counts = engine.get("collector_counts", {})
        positions = counts.get("scanned", 0)
        cycles = counts.get("graph_cycles", 0)
        coverage[row["fixture"]].add(positions)
        cells.append(
            dict(
                fixture=row["fixture"],
                arm=row["arm"],
                complete=row["complete"],
                timing_valid=row["timing_valid"],
                reason=row["reason"],
                positions=positions,
                seconds=row["seconds"],
                work=row["work"],
                counters=counts,
                filter_calls=engine.get("filter_calls", 0),
                pair_dependencies=engine.get("pair_dependencies", 0),
                atoms_per_million_positions=None
                if not positions
                else 10**6 * counts.get("admitted_atoms", 0) / positions,
                cycles_per_million_positions=None
                if not positions
                else 10**6 * cycles / positions,
                atoms_per_cycle=None
                if not cycles
                else counts.get("cycle_atoms", 0) / cycles,
            )
        )
    return dict(
        cells=cells,
        work=4096 * len(results),
        matched_physical_coverage={
            key: len(values) == 1 for key, values in coverage.items()
        },
        interpretation="early cycle proxies only",
    )
