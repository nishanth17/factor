"""Local-only summaries of vmprof samples and streaming PyPy heap records."""

import argparse
import json
from array import array
from collections import Counter, defaultdict
from pathlib import Path


def heap_census(path):
    """Count objects/types without retaining the address/reference graph."""
    totals = defaultdict(lambda: [0, 0, 0])
    header, references, roots, markers, objects = [], 0, 0, 0, 0
    with path.open("rb") as stream:
        while chunk := stream.read(2**20):
            if len(chunk) % 8:
                raise ValueError("truncated heap machine word")
            words = array("q")
            if words.itemsize != 8:
                raise ValueError("this capture requires native 64-bit words")
            words.frombytes(chunk)
            for word in words:
                if len(header) < 3:
                    header.append(word)
                elif word == -1:
                    address, typeindex, size = header
                    if address == 0:
                        if typeindex or size or references or markers:
                            raise ValueError("invalid heap root marker")
                        markers += 1
                    else:
                        if typeindex < 0 or size < 0:
                            raise ValueError("invalid heap object header")
                        row = totals[typeindex]
                        row[0] += 1
                        row[1] += size
                        row[2] += references
                        objects += 1
                        roots += markers == 0
                    header, references = [], 0
                else:
                    references += 1
    if header or markers != 1:
        raise ValueError("truncated heap record or absent root marker")
    names = {}
    typeids = path.parent / "typeids.txt"
    if typeids.exists():
        for line in typeids.read_text().splitlines():
            member, name = line.split(None, 1)
            names[int(member.removeprefix("member"))] = name
    rows = [
        dict(
            typeindex=key,
            name=names.get(key, "unknown"),
            objects=value[0],
            bytes=value[1],
            references=value[2],
        )
        for key, value in totals.items()
    ]
    return dict(
        objects=objects,
        root_objects=roots,
        object_bytes=sum(row[1] for row in totals.values()),
        references=sum(row[2] for row in totals.values()),
        top_types=sorted(rows, key=lambda r: r["bytes"], reverse=True)[:40],
        limitations=(
            "GC object sizes exclude some JIT/native/arena memory. The census "
            "does not compute retained size or identify instance classes."
        ),
    )


def cpu_summary(path):
    import vmprof

    stats = vmprof.read_profile(str(path))
    leaves = Counter()
    unresolved = 0
    for trace, weight, *_ in stats.profiles:
        names = [stats.adr_dict[a] for a in trace if a in stats.adr_dict]
        if names:
            leaves[names[-1]] += weight
        else:
            unresolved += weight
    return dict(
        samples=stats.n_samples,
        expected_samples=stats.expected_samples,
        lost_fraction=stats.get_lost_fraction(),
        inclusive=sorted(
            stats.top_profile(), key=lambda r: r[1], reverse=True
        ),
        last_resolved_function=leaves.most_common(),
        wholly_unresolved_weight=unresolved,
        limitations=(
            "Inclusive functions overlap. Last resolved function can contain "
            "unresolved JIT/native descendants; these are CPU samples, not "
            "wall-time or allocation counts."
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("cpu", "heap"))
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    summary = (
        cpu_summary(args.input)
        if args.kind == "cpu"
        else heap_census(args.input)
    )
    with args.output.open("x") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
