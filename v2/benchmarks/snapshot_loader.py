"""Load the frozen pre-optimization v2 sources in memory for comparisons."""

import hashlib
import json
import sys
import types
from pathlib import Path


def load_snapshot(path=None):
    """Verify a local M8 snapshot and import it as an isolated package.

    Source stays in a provenance JSON artifact; no second checkout or staged
    Python files are created. Call only with the known project-owned snapshot.
    """
    if path is None:
        path = Path(__file__).resolve().parents[1] / "audit"
        path /= "m8_source_snapshot.json"
    data = json.loads(Path(path).read_text())
    name = "_factor_m8"
    package = types.ModuleType(name)
    package.__path__ = []
    package.__package__ = name
    sys.modules[name] = package
    for module_name in (
        "__init__",
        "constants",
        "utils",
        "prime_sieve",
        "pollard_rho",
        "pollard_pm1",
        "ecm",
        "factor",
    ):
        source = data["sources"][module_name]
        digest = hashlib.sha256(source.encode()).hexdigest()
        if digest != data["source_sha256"][module_name]:
            raise ValueError(f"corrupt baseline snapshot: {module_name}")
        if module_name == "__init__":
            module = package
        else:
            qualified = f"{name}.{module_name}"
            module = types.ModuleType(qualified)
            module.__package__ = name
            sys.modules[qualified] = module
            setattr(package, module_name, module)
        # The virtual filename identifies the captured source in tracebacks.
        module.__file__ = str(Path(path).parent / f"M8:{module_name}.py")
        exec(compile(source, module.__file__, "exec"), module.__dict__)
    return package


def load_stage_jobs():
    """Load owned M12 candidates and their frozen arithmetic dependencies.

    Every source byte is hash-checked. Supplied cursor/budget objects remain
    explicit caller controls, allowing committed-boundary comparisons with
    the current implementation without replacing immutable baseline code.
    """
    root = Path(__file__).resolve().parents[1]
    path = root / "audit/m12_source_snapshot.json"
    data = json.loads(path.read_text())
    for name, source in data["sources"].items():
        expected = data["source_sha256"][name]
        if hashlib.sha256(source.encode()).hexdigest() != expected:
            raise ValueError(f"corrupt M12 snapshot: {name}")
    # Freeze the dependencies too. Current additive schedule APIs must not
    # silently alter an immutable candidate control or make it unloadable.
    package_name = "_factor_m12"
    package = types.ModuleType(package_name)
    package.__path__ = []
    package.__package__ = package_name
    sys.modules[package_name] = package
    for name in (
        "constants",
        "utils",
        "prime_sieve",
        "ecm",
        "budget",
        "schedules",
        "stage_jobs",
    ):
        qualified = package_name + "." + name
        module = types.ModuleType(qualified)
        module.__package__ = package_name
        module.__file__ = str(path.parent / f"M12:{name}.py")
        sys.modules[qualified] = module
        setattr(package, name, module)
        exec(
            compile(data["sources"][name], module.__file__, "exec"),
            module.__dict__,
        )
    return package.stage_jobs
