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
    """Load owned M12 candidate code in memory with unchanged dependencies.

    This control isolates batch-loop changes without another checkout. Reject
    modified dependency modules rather than attributing their effects to the
    stage-job implementation. Snapshot hashes identify executable evidence.
    """
    root = Path(__file__).resolve().parents[1]
    path = root / "audit/m12_source_snapshot.json"
    data = json.loads(path.read_text())
    for name, source in data["sources"].items():
        expected = data["source_sha256"][name]
        if hashlib.sha256(source.encode()).hexdigest() != expected:
            raise ValueError(f"corrupt M12 snapshot: {name}")
        if (
            name != "stage_jobs"
            and hashlib.sha256((root / f"{name}.py").read_bytes()).hexdigest()
            != expected
        ):
            raise ValueError(f"M12 dependency changed: {name}")
    module = types.ModuleType("_factor_m12_stage_jobs")
    module.__package__ = "v2"
    module.__file__ = str(path.parent / "M12:stage_jobs.py")
    exec(
        compile(data["sources"]["stage_jobs"], module.__file__, "exec"),
        module.__dict__,
    )
    return module
