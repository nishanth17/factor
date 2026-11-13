"""Common-grant completion quantiles and independent-modulus uncertainty."""

import math
import random

from ..b12_f4.selection import censored_quantile


def interval(clusters, *, grant, draws=4000, seed=2026110423):
    """Keep unfinished early refusals in the full cohort at the common grant.

    Early exhausted searches are informative incompletions, not independent
    censoring times. Removing them from a Kaplan–Meier risk set can invent a
    fast tail. Repeats stay inside their independent modulus cluster.
    """
    if (
        len(clusters) < 12
        or type(draws) is not int
        or not 100 <= draws <= 4000
    ):
        raise ValueError(
            "at least twelve independent moduli and finite draws required"
        )
    if not math.isfinite(grant) or grant <= 0:
        raise ValueError("finite common grant required")
    rows = list(clusters.values())
    for pair in rows:
        if (
            set(pair) != {"control", "nominee"}
            or len(pair["control"]) != 9
            or len(pair["nominee"]) != 9
        ):
            raise ValueError(
                "nine repeats required; stopped cohorts remain inconclusive"
            )
        for observations in pair.values():
            for seconds, complete in observations:
                if (
                    type(complete) is not bool
                    or not math.isfinite(seconds)
                    or not 0 <= seconds <= grant
                ):
                    raise ValueError("invalid common-grant observation")

    def summaries(selected):
        result = {}
        for arm in ("control", "nominee"):
            observations = [
                dict(seconds=seconds, complete=complete, grant=grant)
                for pair in selected
                for seconds, complete in pair[arm]
            ]
            result[arm] = dict(
                median=censored_quantile(observations, 0.5),
                p90=censored_quantile(observations, 0.9),
                completion=sum(row["complete"] for row in observations)
                / len(observations),
            )
        return result

    generator = random.Random(seed)
    samples = {"median": [], "p90": [], "completion": []}
    unidentified = {"median": 0, "p90": 0}
    for _ in range(draws):
        summary = summaries(generator.choices(rows, k=len(rows)))
        for metric in ("median", "p90"):
            a, b = [
                summary[arm][metric]["seconds"]
                for arm in ("control", "nominee")
            ]
            if a is None or b is None:
                unidentified[metric] += 1
            else:
                samples[metric].append(a - b)
        samples["completion"].append(
            summary["nominee"]["completion"] - summary["control"]["completion"]
        )
    intervals = {}
    for metric, values in samples.items():
        values.sort()
        intervals[metric] = (
            None
            if unidentified.get(metric, 0)
            else [
                values[int(0.025 * draws)],
                values[min(draws - 1, int(0.975 * draws))],
            ]
        )
    losses = sum(
        any(a[1] and not b[1] for a, b in zip(row["control"], row["nominee"]))
        for row in rows
    )
    return dict(
        independent_moduli=len(rows),
        estimates=summaries(rows),
        intervals=intervals,
        unidentified_bootstrap_draws=unidentified,
        moduli_with_completion_loss=losses,
        zero_loss_upper_95=1 - 0.05 ** (1 / len(rows))
        if losses == 0
        else None,
        guardrail=(
            "inconclusive unless a reviewed paired zero-loss decision passes"
        ),
    )
