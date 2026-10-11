"""Fresh C9 confirmation with frozen selection and charged paired calls."""

import argparse
import hashlib
import json
import os
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

from ...ecm.p52.p52_realistic import check_quiet
from ...infrastructure.performance.performance_audit import verify_corpus
from ...suites.build_phase_two_corpus import (
    certified_prime,
    verify_certificates,
)
from ..b1.b1_calibration import save
from ..c1.c1_confirmation_support import interval
from . import c9_study as study

SELECTED = study.INPUTS / "controls/c9_selected.json"
CORPUS = study.INPUTS / "corpora/c9_fresh_corpus.json"
ARMS = ("slp", "dlp", "challenger")
SEEDS = (7, 29, 47)
CAPS = {30: 5, 40: 30, 50: 240}
GENERATION_SEED = 202610100909


def hashes():
    result = study.hashes()
    for path in (study.ROOT / "v2/benchmarks").rglob("*.py"):
        if "results" in path.relative_to(study.ROOT).parts:
            continue
        result[str(path.relative_to(study.ROOT))] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
    return result


def load():
    selection = json.loads(SELECTED.read_text())
    if selection["source"] != hashes():
        raise ValueError("C9 confirmation source changed after freeze")
    corpus = json.loads(CORPUS.read_text())
    verify_corpus(corpus)
    if (
        corpus["selection_sha256"]
        != hashlib.sha256(SELECTED.read_bytes()).hexdigest()
    ):
        raise ValueError("corpus does not follow this selection")
    return selection, corpus


class GenerationRandom(random.Random):
    """Bound the certified, intentionally non-RSA sampling population."""

    def __init__(self):
        super().__init__(GENERATION_SEED)
        self.calls = 0
        self.started = time.monotonic()

    def randrange(self, *args, **kwargs):
        self.calls += 1
        if self.calls > 1000000 or time.monotonic() - self.started > 120:
            raise RuntimeError("C9 corpus generation cap")
        return super().randrange(*args, **kwargs)


def existing_inputs():
    """Scan versioned corpus n fields without feeding factors to algorithms."""
    found = set()

    def walk(value):
        if isinstance(value, dict):
            if type(value.get("n")) is int:
                found.add(value["n"])
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)

    for path in (study.INPUTS / "corpora").glob("*.json"):
        walk(json.loads(path.read_text()))
    return found


def generate():
    if CORPUS.exists():
        raise ValueError("C9 fresh corpus already exists; no regeneration")
    selection = json.loads(SELECTED.read_text())
    if selection["source"] != hashes():
        raise ValueError("selection must freeze source before generation")
    seen = existing_inputs()
    original_count = len(seen)
    generator, certificates, fixtures = GenerationRandom(), {}, []
    for digits, count in ((30, 1), (40, 3), (50, 3)):
        bits = (digits * 3322 // 1000 + 1) // 2
        for index in range(count):
            for _ in range(10000):
                p = certified_prime(bits, generator, certificates)
                q = certified_prime(bits, generator, certificates)
                n = p * q
                if p != q and n not in seen and len(str(n)) == digits:
                    break
            else:
                raise RuntimeError("C9 product generation cap")
            seen.add(n)
            fixtures.append(
                dict(
                    id=f"c9_fresh_{digits}_{index}",
                    n=n,
                    factors=sorted((p, q)),
                    kind="balanced",
                    digits=digits,
                    split="held_out",
                )
            )
    verify_certificates(certificates)
    corpus = dict(
        schema=1,
        seed=GENERATION_SEED,
        fixtures=fixtures,
        certificates=certificates,
        generation_calls=generator.calls,
        generation_seconds=time.monotonic() - generator.started,
        disjoint_from_inputs=original_count,
        selection_sha256=hashlib.sha256(SELECTED.read_bytes()).hexdigest(),
        sampling="Balanced certified primes conditioned on a large "
        "p-1 factor; not RSA-distributed. Three independent inputs "
        "per 40/50-digit class; one 30-digit regression input.",
    )
    verify_corpus(corpus)
    save(CORPUS, corpus)


class Expenses:
    """Persist inclusive costs after every call; reserve interrupted calls."""

    def __init__(self, output, selection, corpus):
        self.path = output / "expenses.json"
        if self.path.exists():
            self.data = json.loads(self.path.read_text())
        else:
            self.data = dict(
                prior_wall=selection["spent_before_confirmation"]["wall"],
                prior_cpu=selection["spent_before_confirmation"]["cpu"],
                generation_seconds=corpus["generation_seconds"],
                calls=[],
            )

    def used(self):
        return (
            max(
                self.data["prior_wall"]
                + sum(r["wall"] for r in self.data["calls"]),
                self.data["prior_cpu"]
                + sum(r["cpu"] for r in self.data["calls"]),
            )
            + self.data["generation_seconds"]
        )

    def write(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".tmp")
        temporary.write_text(json.dumps(self.data, indent=2) + "\n")
        temporary.replace(self.path)

    def call(self, label, cap, function):
        if self.used() + cap > 14400:
            raise RuntimeError("C9 aggregate allowance cannot cover next call")
        if (
            sum(
                p.stat().st_size
                for p in self.path.parent.rglob("*")
                if p.is_file()
            )
            > 256 * 2**20
        ):
            raise RuntimeError("C9 capture disk cap")
        record = dict(label=label, wall=cap, cpu=cap, status="reserved")
        self.data["calls"].append(record)
        self.write()
        start, cpu = time.monotonic(), time.process_time()
        child = resource.getrusage(resource.RUSAGE_CHILDREN)
        child_before = child.ru_utime + child.ru_stime
        try:
            result = function()
            record["status"] = "complete"
            return result
        except BaseException:
            record["status"] = "failed"
            raise
        finally:
            child = resource.getrusage(resource.RUSAGE_CHILDREN)
            record.update(
                wall=time.monotonic() - start,
                cpu=time.process_time()
                - cpu
                + child.ru_utime
                + child.ru_stime
                - child_before,
            )
            self.write()


def cold(fixture, arm, output):
    path = output / f"cold-{fixture['digits']}-{arm}.json"
    start = time.monotonic()
    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            "-m",
            __spec__.name,
            "cold-worker",
            "--output",
            str(path),
            "--fixture",
            fixture["id"],
            "--arm",
            arm,
        ],
        capture_output=True,
        text=True,
        timeout=CAPS[fixture["digits"]] + 30,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr)
    save(
        path.with_name(path.stem + "-startup.json"),
        dict(
            process_seconds=time.monotonic() - start,
            stdout=completed.stdout,
            stderr=completed.stderr,
            classification="One cold diagnostic, separate from warm samples",
        ),
    )


def unstable(rows):
    """Extend conflicting paired signs or repeated drift, as prespecified."""
    pairs = {}
    for row in rows:
        pairs.setdefault((row["id"], row["sample"]), {})[row["arm"]] = row
    reasons = []
    for reference in ("slp", "dlp"):
        for identity in {key[0] for key in pairs}:
            effects = [
                1
                - pair["challenger"]["capped_seconds"]
                / pair[reference]["capped_seconds"]
                for (case, _), pair in pairs.items()
                if case == identity and set(pair) == set(ARMS)
            ]
            if effects and min(effects) < 0 < max(effects):
                reasons.append(
                    "within-input paired direction differs versus "
                    + reference
                    + ": "
                    + identity
                )
    cells = {}
    for row in rows:
        cells.setdefault((row["id"], row["seed"], row["arm"]), []).append(row)
    for key, values in cells.items():
        if len(values) < 2:
            continue
        values = [
            r["seconds"] for r in sorted(values, key=lambda r: r["sample"])
        ]
        median = statistics.median(values)
        if abs(values[-1] - values[0]) / median > 0.10:
            reasons.append("repeated-cell drift: " + str(key))
        if len(values) >= 4:
            q = statistics.quantiles(values, n=4, method="inclusive")
            if (q[2] - q[0]) / median > 0.10:
                reasons.append("repeated-cell IQR: " + str(key))
    return reasons


def confirm(output, digits):
    selection, corpus = load()
    metadata = output / "manifest.json"
    expected = dict(
        source=selection["source"],
        runtime=sys.version,
        platform=platform.platform(),
        selected_sha256=hashlib.sha256(SELECTED.read_bytes()).hexdigest(),
        corpus_sha256=hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
    )
    if metadata.exists():
        if json.loads(metadata.read_text()) != expected:
            raise ValueError(
                "confirmation output belongs to a different freeze"
            )
    else:
        save(metadata, expected)
    if (output / f"finished-{digits}.json").exists():
        # Repeating a finished band cannot acquire additional samples.
        return
    expense = Expenses(output, selection, corpus)
    fixtures = [f for f in corpus["fixtures"] if f["digits"] == digits]
    configs = {
        a: study.decode(selection["configurations"][str(digits)][a])
        for a in ARMS
    }
    for arm in ARMS:
        cold_path = output / f"cold-{digits}-{arm}.json"
        if not cold_path.exists():
            expense.call(
                f"cold-{digits}-{arm}",
                CAPS[digits] + 30,
                lambda a=arm: cold(fixtures[0], a, output),
            )
        path = output / f"warmup-{digits}-{arm}.json"
        if not path.exists():
            warm, elapsed, elapsed_cpu = [], 0.0, 0.0
            while min(elapsed, elapsed_cpu) < 3:
                row = expense.call(
                    f"warmup-{digits}-{arm}",
                    5,
                    lambda a=arm: study.run_one(
                        fixtures[0], 7, configs[a], a, 3.1
                    ),
                )
                warm.append(row)
                elapsed += row["seconds"]
                elapsed_cpu += row["cpu_seconds"]
            save(path, warm)
    rows = []
    for target in (9, 18):
        for index in range(len(rows) // len(ARMS), target):
            fixture = fixtures[index % len(fixtures)]
            sample = index // len(fixtures)
            seed = SEEDS[sample % 3] if digits != 30 else SEEDS[index % 3]
            order = list(ARMS)
            random.Random(202610109 + digits * 100 + index).shuffle(order)
            for arm in order:
                path = output / f"{fixture['id']}-{sample:02d}-{arm}.json"
                if path.exists():
                    row = json.loads(path.read_text())
                else:
                    row = expense.call(
                        path.stem,
                        CAPS[digits] + 2,
                        lambda: study.run_one(
                            fixture, seed, configs[arm], arm, CAPS[digits]
                        ),
                    )
                    row["sample"] = sample
                    save(path, row)
                    print(
                        json.dumps(
                            {
                                k: row[k]
                                for k in (
                                    "id",
                                    "arm",
                                    "sample",
                                    "complete",
                                    "seconds",
                                )
                            }
                        ),
                        flush=True,
                    )
                rows.append(row)
        reasons = unstable(rows)
        # Identical smaller-band DLP/challenger configurations are regression
        # controls, not a separate speed claim. Do not extend random signs.
        if configs["challenger"] == configs["dlp"]:
            reasons = [r for r in reasons if r.startswith("repeated-cell")]
        save(
            output / f"stability-{digits}-{target}.json",
            dict(
                reasons=reasons,
                matched_samples=target,
                further_extension=target == 9 and bool(reasons),
            ),
        )
        if not reasons:
            break
        cached_extension = sum(
            1
            for path in output.glob("c9_fresh_*.json")
            if (saved := json.loads(path.read_text()))["digits"] == digits
            and saved["sample"] >= 9 // len(fixtures)
        )
        extension_cap = max(0, 27 - cached_extension) * (CAPS[digits] + 2)
        if target == 9 and expense.used() + extension_cap > 14400:
            save(
                output / f"extension-refused-{digits}.json",
                dict(
                    reason="aggregate cap",
                    used=expense.used(),
                    reasons=reasons,
                    verdict="timing inconclusive; retain baseline",
                ),
            )
            break
    save(
        output / f"finished-{digits}.json",
        dict(
            matched_samples=len(rows) // len(ARMS),
            total_charged=expense.used(),
        ),
    )


def arm_summary(rows):
    """Keep resource, collection and dependency outcomes beside cost ratios."""

    def median(values):
        return statistics.median(values)

    engines = [r["stats"].get("engine", {}) for r in rows]
    counts = [e.get("collector_counts", {}) for e in engines]
    filters = [e.get("last_filter", {}) for e in engines]
    return dict(
        calls=len(rows),
        complete=sum(r["complete"] for r in rows),
        proper_factors=sum(r["divisor"] is not None for r in rows),
        wall_total=sum(r["seconds"] for r in rows),
        cpu_total=sum(r["cpu_seconds"] for r in rows),
        capped_total=sum(r["capped_seconds"] for r in rows),
        wall_median=median(r["seconds"] for r in rows),
        work_total=sum(r["work"] for r in rows),
        reasons=sorted({r["reason"] for r in rows}),
        peak_owned_bytes=max(
            r["stats"].get("workspace_bytes", 0) for r in rows
        ),
        process_peak_rss=max(r["rss_bytes"] for r in rows),
        median_counters={
            key: median(c.get(key, 0) for c in counts)
            for key in (
                "scanned",
                "candidates",
                "refined_rejections",
                "admitted_atoms",
                "duplicates",
                "matches",
                "evictions",
                "split_attempts",
                "rho_evaluations",
                "split_failed",
                "endpoint_bound",
                "product_rejections",
                "prime_rejections",
                "not_two_primes",
                "split_quota_rejections",
                "graph_cycles",
                "dlp_cycles",
                "cycle_atoms",
                "long_cycle_rejections",
                "graph_rebuilds",
            )
        },
        median_filter={
            key: median(f.get(key, 0) for f in filters)
            for key in (
                "input_rows",
                "output_rows",
                "output_columns",
                "singletons_removed",
                "weight_two_merges",
                "zero_dependencies",
            )
        },
        total_stage_seconds={
            stage: sum(
                e.get("stage_seconds", {}).get(stage, 0) for e in engines
            )
            for stage in (
                "collection",
                "preparation",
                "filtering",
                "elimination",
                "extraction",
            )
        },
        last_dependency_counts=[
            e.get("last_solver", {}).get("dependencies", 0) for e in engines
        ],
        trivial_dependencies=sum(
            e.get("trivial_dependencies", 0) for e in engines
        ),
        filter_calls=sum(e.get("filter_calls", 0) for e in engines),
        solve_calls=sum(e.get("solve_calls", 0) for e in engines),
    )


def summarize(output):
    groups = {}
    for path in sorted(output.glob("c9_fresh_*.json")):
        row = json.loads(path.read_text())
        groups.setdefault(row["digits"], {}).setdefault(
            row["id"], {}
        ).setdefault(row["sample"], {})[row["arm"]] = row
    result = {}
    for digits, cases in groups.items():
        all_rows = [
            r
            for samples in cases.values()
            for pair in samples.values()
            for r in pair.values()
        ]
        band = {
            "arms": {
                arm: arm_summary([r for r in all_rows if r["arm"] == arm])
                for arm in ARMS
            }
        }
        for reference, candidate in (
            ("slp", "challenger"),
            ("dlp", "challenger"),
            ("dlp", "slp"),
        ):
            paired = [
                [
                    (p[reference], p[candidate])
                    for p in samples.values()
                    if set(p) == set(ARMS)
                ]
                for samples in cases.values()
            ]
            paired = [p for p in paired if p]
            generator, benefits = random.Random(90900 + digits), []
            for _ in range(10000):
                old = new = 0
                for _ in paired:
                    source = generator.choice(paired)
                    for _ in source:
                        a, b = generator.choice(source)
                        old += a["capped_seconds"]
                        new += b["capped_seconds"]
                benefits.append(1 - new / old)
            flat = [pair for source in paired for pair in source]
            band[candidate + "_vs_" + reference] = dict(
                reference=reference,
                candidate=candidate,
                benefit=1
                - sum(b["capped_seconds"] for a, b in flat)
                / sum(a["capped_seconds"] for a, b in flat),
                interval95=interval(benefits),
                independent_inputs=len(paired),
                matched_samples=len(flat),
                complete_reference=sum(a["complete"] for a, b in flat),
                complete_candidate=sum(b["complete"] for a, b in flat),
                total_reference=sum(a["capped_seconds"] for a, b in flat),
                total_candidate=sum(b["capped_seconds"] for a, b in flat),
                per_input={
                    pairs[0][0]["id"]: dict(
                        benefit=1
                        - sum(b["capped_seconds"] for a, b in pairs)
                        / sum(a["capped_seconds"] for a, b in pairs),
                        samples=len(pairs),
                        complete_reference=sum(
                            a["complete"] for a, b in pairs
                        ),
                        complete_candidate=sum(
                            b["complete"] for a, b in pairs
                        ),
                    )
                    for pairs in paired
                },
            )
        result[str(digits)] = band
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode", choices=("generate", "confirm", "cold-worker", "summarize")
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--digits", type=int, choices=(30, 40, 50))
    parser.add_argument("--fixture")
    parser.add_argument("--arm", choices=ARMS)
    args = parser.parse_args()
    study.runtime()
    if args.mode == "cold-worker":
        selection, corpus = load()
        fixture = next(
            f for f in corpus["fixtures"] if f["id"] == args.fixture
        )
        config = study.decode(
            selection["configurations"][str(fixture["digits"])][args.arm]
        )
        with patch.object(
            study,
            "check_quiet",
            lambda allowed: check_quiet(allowed | {os.getppid()}),
        ):
            save(
                args.output,
                study.run_one(
                    fixture, 7, config, args.arm, CAPS[fixture["digits"]]
                ),
            )
    else:
        with study.machine_window():
            if args.mode == "generate":
                generate()
            elif args.mode == "confirm":
                if args.output is None or args.digits is None:
                    parser.error("confirmation requires output and digits")
                confirm(args.output, args.digits)
            else:
                save(args.output / "summary.json", summarize(args.output))


if __name__ == "__main__":
    main()
