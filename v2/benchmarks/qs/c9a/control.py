"""Load the immutable committed control without relying on another worktree."""

import hashlib
import importlib
import importlib.util
import json
import sys
import tempfile
import zipfile
from pathlib import Path

INPUTS = Path(__file__).resolve().parents[2] / "inputs"
ARCHIVE = INPUTS / "controls/c9a_mainline.zip"
MANIFEST = INPUTS / "controls/c9a_mainline.json"
_DIRECTORY = None


def control_modules():
    """Verify all bytes before importing a private, process-local package."""
    global _DIRECTORY
    name = "c9a_immutable_control"
    if name not in sys.modules:
        manifest = json.loads(MANIFEST.read_text())
        if (
            hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
            != (manifest["archive_sha256"])
        ):
            raise ValueError("immutable control archive changed")
        _DIRECTORY = tempfile.TemporaryDirectory(prefix="factor-c9a-control-")
        root = Path(_DIRECTORY.name)
        with zipfile.ZipFile(ARCHIVE) as archive:
            if set(archive.namelist()) != set(manifest["files"]):
                raise ValueError("immutable control inventory changed")
            for path, expected in manifest["files"].items():
                parts = Path(path).parts
                if Path(path).is_absolute() or ".." in parts:
                    raise ValueError("unsafe control archive path")
                data = archive.read(path)
                if hashlib.sha256(data).hexdigest() != expected:
                    raise ValueError("immutable control member changed")
                target = root / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
        package = root / "v2"
        spec = importlib.util.spec_from_file_location(
            name,
            package / "__init__.py",
            submodule_search_locations=[str(package)],
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return (
        importlib.import_module(name + ".qs"),
        importlib.import_module(name + ".execution.budget"),
    )
