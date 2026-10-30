"""Strict owned-source replay of committed selected A1a integration."""

import hashlib
import importlib
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASELINE = ROOT / "v2/benchmarks/inputs/baselines/a1a_extension_mainline.json"
START = "fee6a41c9af2702c34f964536e71ad3011d4cb22"
PACKAGE = "_a1a_extension_control"


def checked_baseline():
    """Reject changed source text, paths and baseline provenance."""
    baseline = json.loads(BASELINE.read_text())
    if baseline["schema"] != 1 or baseline["commit"] != START:
        raise ValueError("unexpected prior A1a source identity")
    if baseline["source"].keys() != baseline["source_sha256"].keys():
        raise ValueError("incomplete prior A1a source hashes")
    for name, source in baseline["source"].items():
        path = Path(name)
        if (
            path.is_absolute()
            or path.parts[0] != "v2"
            or any(part in (".", "..") for part in path.parts)
            or path.suffix not in (".py", ".json")
            or hashlib.sha256(source.encode()).hexdigest()
            != baseline["source_sha256"][name]
        ):
            raise ValueError("invalid or corrupt owned baseline source")
    return baseline


def load_control(output):
    """Load control classes separately; never share live relation types."""
    baseline = checked_baseline()
    directory = output / "control_runtime" / PACKAGE
    for name, source in baseline["source"].items():
        destination = directory.joinpath(*Path(name).parts[1:])
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if destination.read_text() != source:
                raise ValueError("materialized control source changed")
        else:
            destination.write_text(source)

    if PACKAGE in sys.modules:
        if Path(sys.modules[PACKAGE].__file__).parent != directory:
            raise ValueError("control package already belongs to another run")
    else:
        spec = importlib.util.spec_from_file_location(
            PACKAGE,
            directory / "__init__.py",
            submodule_search_locations=[str(directory)],
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[PACKAGE] = module
        spec.loader.exec_module(module)
    return (
        importlib.import_module(PACKAGE + ".qs"),
        importlib.import_module(PACKAGE + ".execution.budget"),
        importlib.import_module(PACKAGE + ".common.utils"),
    )
