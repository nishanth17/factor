"""Import owned pre-P3.3 sources in memory for matched causal comparisons."""

import hashlib
import json
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "audit/m25_p33_before_sources.json"


def load_qs_arm(name, changes=()):
    """Load old QS plus selected current modules and common postprocessing.

    Sources remain directly readable provenance, not a staging checkout.
    Every old source hash is checked before execution. New postprocessing
    uses each arm's own relation classes, avoiding type/conversion asymmetry.
    """
    data = json.loads(SNAPSHOT.read_text())
    sources = data["sources"]
    for module_name, source in sources.items():
        if (
            hashlib.sha256(source.encode()).hexdigest()
            != (data["source_sha256"][module_name])
        ):
            raise ValueError("corrupt owned QS source snapshot")
    package = types.ModuleType(name)
    package.__path__ = []
    package.__package__ = name
    sys.modules[name] = package
    qs = types.ModuleType(name + ".qs")
    qs.__path__ = []
    qs.__package__ = name + ".qs"
    sys.modules[qs.__name__] = qs
    package.qs = qs
    names = [key for key in sources if key != "qs.__init__"]
    names += ["qs.linear_algebra", "qs.extraction", "qs.pipeline"]
    for module_name in names:
        current = module_name in changes or module_name not in sources
        path = ROOT.joinpath(*module_name.split(".")).with_suffix(".py")
        source = path.read_text() if current else sources[module_name]
        if module_name == "__init__":
            module = package
        else:
            qualified = name + "." + module_name
            module = types.ModuleType(qualified)
            module.__package__ = qualified.rsplit(".", 1)[0]
            sys.modules[qualified] = module
            parent_name, _, child_name = qualified.rpartition(".")
            setattr(sys.modules[parent_name], child_name, module)
        module.__file__ = (
            str(path if current else SNAPSHOT) + ":" + module_name
        )
        exec(compile(source, module.__file__, "exec"), module.__dict__)
    return package
