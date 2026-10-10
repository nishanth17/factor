"""Supplement R3 with mixed full/matched payloads and live coexistence caps."""

import argparse
import hashlib
import importlib
import json
import platform
import sys
from pathlib import Path

from ...infrastructure.performance.performance_audit import (
    fingerprint,
    verify_corpus,
)
from ...support.paths import (
    source_path,
)
from .p38_r3 import (
    HERE,
    MEMORY,
    ROOT,
    WORK,
    allowance,
    history_filter,
    live_compaction,
    modules,
    paired_measure,
    prepare_fixture,
    rank_oracle,
    rebuilt_filter,
)


def run(mods, prepared, variant, rank):
    """Include incidence, fill, lifting and every exact modular extraction."""
    algebra = mods["qs.linear_algebra"]
    budget = allowance(mods["budget"])
    remaining = MEMORY - prepared.workspace_bytes
    if variant == "rebuild":
        matrix = rebuilt_filter(algebra)(
            prepared.rows,
            weight_two=True,
            budget=budget,
            memory_bytes=remaining,
        )
    elif variant.startswith(("history", "dense_batch")):
        matrix = history_filter(
            prepared.rows,
            algebra,
            budget,
            remaining,
            batch_size=32 if variant.endswith("32") else 1,
            use_history=variant.startswith("history"),
        )
    else:
        matrix = algebra.filter_matrix(
            prepared.rows,
            weight_two=True,
            budget=budget,
            memory_bytes=remaining,
        )
        if variant == "live":
            matrix = live_compaction(matrix, budget, remaining)

    dependencies = algebra.DependencySolver(matrix, budget=budget).run()

    if len(dependencies) != len(prepared.rows) - rank:
        raise AssertionError("mixed kernel dimension differs from oracle")
    congruences = []

    for mask in dependencies:
        algebra.verify_dependency(mask, prepared.rows)
        congruences.append(
            mods["qs.extraction"].extract_dependency(
                prepared,
                mask,
                budget=budget,
            )
        )

    owned = prepared.workspace_bytes + matrix.workspace_bytes
    if owned > MEMORY or budget.used > WORK:
        raise AssertionError("mixed provenance/matrix allowance exceeded")
    return dict(
        dependencies=len(dependencies),
        extracted=len(congruences),
        proper_divisors=sum(c.divisor is not None for c in congruences),
        work=budget.used,
        owned_bytes=owned,
        matrix_stats=matrix.stats,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("preserve captures; choose a new output path")
    corpus = json.loads(
        (
            source_path(HERE / "inputs/corpora/p38_r3_training_corpus.json")
        ).read_text()
    )
    verify_corpus(corpus)
    mods = modules(importlib.import_module("v2"))
    before = fingerprint(ROOT)
    results = {}

    for fixture in corpus["fixtures"][:2]:
        base, relations, atoms = prepare_fixture(mods, fixture)
        from ....qs.relations import AtomicRelation, CombinedRelation

        full = [r for r in relations if isinstance(r, AtomicRelation)][:64]
        combined = [r for r in relations if isinstance(r, CombinedRelation)][
            :64
        ]
        if len(full) != 64 or len(combined) != 64:
            raise AssertionError("mixed trace requires 64 full and 64 matched")
        mixed = tuple(r for pair in zip(full, combined) for r in pair)
        prepared = mods["qs.extraction"].prepare_relations(
            mixed,
            base,
            atoms,
            budget=allowance(mods["budget"]),
            memory_bytes=MEMORY,
        )
        rank = rank_oracle(prepared.rows)
        variants = (
            "rebuild",
            "current",
            "live",
            "dense_batch1",
            "dense_batch32",
            "history1",
            "history32",
        )
        results[fixture["id"]] = paired_measure(
            {
                name: lambda name=name: run(mods, prepared, name, rank)
                for name in variants
            }
        )

    if before != fingerprint(ROOT):
        raise AssertionError("runtime changed during mixed capture")
    args.output.write_text(
        json.dumps(
            dict(
                source_sha256=before,
                driver_sha256=hashlib.sha256(
                    source_path(Path(__file__)).read_bytes()
                ).hexdigest(),
                python=platform.python_version(),
                implementation="PyPy",
                trace="64 full and 64 matched rows interleaved; "
                "every mask extracted",
                budgets=dict(work=WORK, wall_cpu=30, owned_bytes=MEMORY),
                results=results,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
