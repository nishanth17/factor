"""Owned current-mainline control for the v4 experiment design."""

import hashlib
import importlib
import importlib.util
import json
import sys
from dataclasses import asdict
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)

from .configuration_v4 import configuration

ROOT = Path(__file__).resolve().parents[4]
BASELINE = ROOT / "v2/benchmarks/inputs/baselines/c10a01_mainline_0076e25.json"
START = "0076e256871722d52d849ab344152e0a7535c2ce"
PACKAGE = "_c10a01_mainline_0076e25"
PARENT_OWNED = 32 * 2**20
PARENT_RSS = 256 * 2**20


def checked_baseline():
    """Check source identity and finite safe paths before materialization."""
    if BASELINE.stat().st_size > 4 * 2**20:
        raise ValueError("owned control exceeds its source allowance")
    baseline = json.loads(BASELINE.read_text())
    if baseline["schema"] != 1 or baseline["commit"] != START:
        raise ValueError("unexpected mainline control identity")
    if baseline["source"].keys() != baseline["source_sha256"].keys():
        raise ValueError("incomplete control source hashes")
    for name, source in baseline["source"].items():
        path = Path(name)
        if (
            path.is_absolute()
            or path.parts[0] != "v2"
            or any(part in (".", "..") for part in path.parts)
            or path.suffix != ".py"
            or hashlib.sha256(source.encode()).hexdigest()
            != baseline["source_sha256"][name]
        ):
            raise ValueError("invalid owned control source")
    return baseline


def load_control(output):
    """Load isolated classes without importing challenger arithmetic."""
    baseline = checked_baseline()
    directory = Path(output) / "control_runtime" / PACKAGE
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
            raise ValueError("control package belongs to another directory")
    else:
        spec = importlib.util.spec_from_file_location(
            PACKAGE,
            directory / "__init__.py",
            submodule_search_locations=[str(directory)],
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[PACKAGE] = module
        spec.loader.exec_module(module)
    return tuple(
        importlib.import_module(PACKAGE + suffix)
        for suffix in (".qs", ".execution.budget", ".common.utils")
    )


def configuration_values(digits):
    """Translate only disabled new fields; preserve every numerical setting."""
    values = asdict(configuration(digits, "control"))
    for name in (
        "polynomial_policy",
        "multiplier_policy",
        "root_update_policy",
    ):
        if values.pop(name) != "standard":
            raise ValueError("challenger policy cannot enter mainline control")
    values["memory_bytes"] -= PARENT_OWNED
    return values


def control_configuration(digits, module):
    """Instantiate the exact literal configuration in the owned runtime."""
    values = configuration_values(digits)
    values["collector"] = adapt_resource_configuration(
        module.GrowingDoubleLargeSieveConfig, **values["collector"]
    )
    return adapt_resource_configuration(module.CountedSIQSConfig, **values)
