"""Prespecified nearest/flyer and legacy controls at the feasible band."""

import argparse
import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

from .build_phase_two_corpus import verify_certificates
from .p38_r1_capacity import decode_config, measure, source_hashes

POLICIES = ("nearest", "flyer", "reference")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        parser.error("PyPy implementing Python 3.11 is required")
    if args.output.exists():
        parser.error("preserve existing results")
    corpus = json.loads(args.corpus.read_text())
    frozen_bytes = args.frozen.read_bytes()
    frozen = json.loads(frozen_bytes)
    digest = hashlib.sha256(frozen_bytes).hexdigest()
    hashes = source_hashes()
    if corpus["frozen_sha256"] != digest or hashes != frozen["source_sha256"]:
        parser.error("confirmation must match the frozen control and runtime")
    if corpus["seeds"] != frozen["seeds"]:
        parser.error("seeds differ from frozen control")
    verify_certificates(corpus["certificates"])
    fixtures = [
        f
        for f in corpus["fixtures"]
        if f["kind"] == "balanced" and f["digits"] == 30
    ]
    config = decode_config(frozen["configs"]["30"])
    results = dict(
        source_sha256=hashes,
        frozen_sha256=digest,
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        corpus_sha256=hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
        scope="Three prespecified A policies at the trained settings. "
        "Reference retains 64 families and full Gray reuse. Nearest/flyer "
        "share streamed quotas. All arms share resource caps.",
        policies={},
    )

    for policy in POLICIES:
        if policy == "reference":
            # Construct directly: a reference config cannot temporarily carry
            # the streaming-only quota, even as an intermediate dataclass.
            arm = replace(
                config,
                assignment_policy=policy,
                family_count=64,
                polynomials_per_family=0,
            )
        else:
            arm = replace(config, assignment_policy=policy)

        results["policies"][policy] = measure(
            fixtures,
            arm,
            "siqs",
            frozen["confirmation_seconds"]["30"],
            corpus["seeds"],
        )
        if hashes != source_hashes():
            raise RuntimeError("source changed during capture")
        args.output.write_text(json.dumps(results, indent=2) + "\n")
        print(policy, "finished", flush=True)


if __name__ == "__main__":
    main()
