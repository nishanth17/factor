"""Package unchanged measured costs with new source pins; no calibration."""

import argparse
import ast
import hashlib
import json
import re
from pathlib import Path

SOURCE_ROOT = Path(__file__).resolve().parents[3]
MAX_PROFILE_BYTES = 256 * 2**10
MAX_SOURCE_BYTES = 2 * 2**20


def build_profile(reference, expected_sha256, *, source_root=SOURCE_ROOT):
    """Transfer estimates explicitly, preserving the immutable reference."""
    with Path(reference).open("rb") as stream:
        raw = stream.read(MAX_PROFILE_BYTES + 1)
    if len(raw) > MAX_PROFILE_BYTES or hashlib.sha256(raw).hexdigest() != (
        expected_sha256
    ):
        raise ValueError("measured reference hash or size differs")
    data = json.loads(raw)
    tree = ast.parse((source_root / "execution/ecm_selection.py").read_text())
    declaration = next(
        node.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "SOURCE_NAMES"
            for target in node.targets
        )
    )
    names = ast.literal_eval(declaration.args[0])
    pins, total = {}, 0
    for name in names:
        path = source_root / name
        total += path.stat().st_size
        if total > MAX_SOURCE_BYTES:
            raise ValueError("integration source exceeds profile ownership")
        pins[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    data["source_sha256"] = pins
    # Repinning code cannot turn historical estimates into new confirmation.
    data["status"] = "training"
    return data


def repin_adapter(text, sha256):
    """Change exactly the packaged-policy constant, never curve identities."""
    if re.fullmatch(r"[0-9a-f]{64}", sha256) is None:
        raise ValueError("invalid policy hash")
    result, count = re.subn(
        r'(POLICY_SHA256 = \(\n    ")[0-9a-f]{64}("\n\))',
        lambda match: match[1] + sha256 + match[2],
        text,
    )
    if count != 1:
        raise ValueError("expected exactly one policy constant")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--adapter", required=True, type=Path)
    args = parser.parse_args()
    if args.reference.resolve() in (
        args.output.resolve(),
        args.adapter.resolve(),
    ):
        raise ValueError("measured reference must remain immutable")
    data = build_profile(args.reference, args.reference_sha256)
    raw = (json.dumps(data, indent=2) + "\n").encode()
    sha256 = hashlib.sha256(raw).hexdigest()
    adapter = repin_adapter(args.adapter.read_text(), sha256)
    args.output.write_bytes(raw)
    args.adapter.write_text(adapter)
    print(sha256)


if __name__ == "__main__":
    main()
