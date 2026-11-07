"""Reserve a diagnostic companion for every measured configuration/input."""

from copy import deepcopy


def paired_cells(cells):
    """Keep clean evidence and whole-execution CPU/live-heap diagnostics.

    The diagnostic repeats one sample with exactly the same configuration,
    input, seed and lifetime. CPU traces cover every owned worker. Parent and
    one worker have live heaps, with explicitly labeled teardown fallback for
    very short calls. Neither diagnostic timing nor heap pauses rank arms.
    """
    result = []
    for ordinal, original in enumerate(cells):
        if original.get("mode", "clean") != "clean":
            raise ValueError(
                "only clean calls can request a diagnostic companion"
            )
        clean = deepcopy(original)
        name = clean.get("variant", clean["arm"])
        pair = f"{clean['fixture']}:{name}:{ordinal}"
        clean.update(mode="clean", profile_required=True, pair=pair)
        diagnostic = deepcopy(clean)
        diagnostic.pop("profile_required")
        diagnostic.update(
            mode="profile",
            full_profile=True,
            diagnostic_for=pair,
            samples=1,
            observe=True,
        )
        result.extend((clean, diagnostic))
    return result
