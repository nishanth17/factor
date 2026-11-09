"""Freeze the finite G2 breadth proposals; construction runs no factorer."""

from .combined import cell, corpus_identity
from .regime_search import candidate
from .runner import source
from .stage_driver import ROWS

WALL = {40: 30, 50: 90, 55: 180, 60: 300, 65: 900, 70: 1800}
COUNTS = {40: 4, 50: 4, 55: 4, 60: 8, 65: 8, 70: 10}


def prepare(identity, reference_source, folder):
    """Retained incumbents are proposals, not new per-regime winners."""
    commit = source()
    corpus = corpus_identity("regimes")
    batches = []
    for digits, lifetime in sorted(ROWS):
        for role in ("mainline", "serial", "process"):
            name = f"{identity}-{digits}-{lifetime}-{role}"
            spec = cell(
                "S0",
                f"regime_calibration_{digits}_0",
                stage="G2",
                digits=digits,
                lifetime=lifetime,
                endpoint="complete",
                samples=1,
                mode="clean",
                observe=True,
                tier_a=False,
                seed=104729,
                wall=WALL[digits],
                cpu=10 * WALL[digits],
                memory=22 * 2**30,
                work=96 * 10**12,
                checkpoint_bytes=64 * 2**20,
                qs_memory_bytes=4 * 2**30,
            )
            if role == "mainline":
                spec.update(arm="MAINLINE", policy="G0")
            elif role == "process":
                spec.update(nursery_mb=4, max_batch_atoms=2048)
                recipe = candidate(
                    COUNTS[digits],
                    16,
                    features=dict(
                        worker_memory_bytes=512 * 2**20,
                        state_cache_bytes=64 * 2**20,
                        packet_bytes=4 * 2**20,
                        coefficient_cache=True,
                    ),
                )
                recipe.pop("name")
                spec.update(recipe)
            batch = dict(
                identity=name,
                source=commit,
                corpus_kind="regimes",
                corpus_sha256=corpus,
                cells=[spec],
                wall=WALL[digits] + 50,
                cpu=10 * WALL[digits] + 60,
                storage=8 * 2**20,
                profile_folder=str(folder / name / "profiles"),
            )
            if role == "mainline":
                batch["reference_source"] = reference_source
            batches.append(batch)
    return dict(
        identity=identity,
        stage="G2",
        source=commit,
        corpus_kind="regimes",
        corpus_sha256=corpus,
        wall=5 * 3600,
        cpu=50 * 3600,
        storage=4 * 2**30,
        recipe_caps={f"{d}:{life}": 3 for d, life in ROWS},
        batches=batches,
    )
