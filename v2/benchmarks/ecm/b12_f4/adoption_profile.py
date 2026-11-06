"""Package unchanged measured costs with new source pins; no calibration."""

import argparse
import hashlib
import json
import re
from pathlib import Path

MAX_PROFILE_BYTES = 256 * 2**10


def build_profile(reference, expected_sha256):
    """Transfer measured costs without turning them into fresh confirmation."""
    from ....ecm.selection import COMPATIBILITY, PROFILE_VERSION

    raw = Path(reference).read_bytes()
    if (
        len(raw) > MAX_PROFILE_BYTES
        or hashlib.sha256(raw).hexdigest() != expected_sha256
    ):
        raise ValueError("measured reference hash or size differs")
    data = json.loads(raw)
    data.pop("source_sha256", None)
    data.update(
        version=PROFILE_VERSION, compatibility=COMPATIBILITY, status="training"
    )
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
