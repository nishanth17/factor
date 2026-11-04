"""Finite heap/profile summaries; no factoring or worker imports."""

import array
import collections
import importlib.util
from pathlib import Path


def heap_summary(path, maximum=32 * 2**20):
    """Check native-word RPython records and summarize GC-managed live types.

    This format is documented by gc.dump_rpy_heap. Raw allocations, allocator
    freelists and JIT code require the separate process/GC statistics.
    """
    path = Path(path)
    size = path.stat().st_size
    values = array.array("l")
    if not size or size > maximum or size % values.itemsize:
        raise ValueError("heap dump is empty, oversized or word-truncated")
    with path.open("rb") as stream:
        values.fromfile(stream, size // values.itemsize)
    counts, sizes = collections.Counter(), collections.Counter()
    live, references = 0, 0
    cursor, markers, objects = 0, 0, 0
    largest = []
    while cursor < len(values):
        if cursor + 3 >= len(values):
            raise ValueError("truncated heap record")
        address = values[cursor]
        kind, owned = values[cursor + 1], values[cursor + 2]
        if address < 0 or kind < 0 or owned < 0:
            raise ValueError("invalid heap object header")
        end = cursor + 3
        while end < len(values) and values[end] != -1:
            end += 1
        if end == len(values):
            raise ValueError("unterminated heap references")
        if address == 0:
            if (kind, owned, end - cursor) != (0, 0, 3):
                raise ValueError("invalid roots delimiter")
            markers += 1
        else:
            objects += 1
            counts[kind] += 1
            sizes[kind] += owned
            live += owned
            references += end - cursor - 3
            if owned >= 65536:
                largest.append((owned, kind))
        cursor = end + 1
    if markers != 1:
        raise ValueError("missing or duplicate roots delimiter")
    names = {}
    typeids = path.parent / "typeids.txt"
    if typeids.exists():
        names = dict(enumerate(typeids.read_text().splitlines()))
    return dict(
        file=str(path),
        bytes=size,
        objects=objects,
        references=references,
        gc_object_bytes=live,
        types=[
            dict(
                typeid=k,
                name=names.get(k, "unknown"),
                count=counts[k],
                bytes=value,
            )
            for k, value in sizes.most_common(25)
        ],
        largest=[
            dict(bytes=value, typeid=k, name=names.get(k, "unknown"))
            for value, k in sorted(largest, reverse=True)[:20]
        ],
        coverage="GC-managed records; not total process RSS or all raw memory",
    )


def cpu_summary(path, reader_path, maximum=32 * 2**20):
    """Use a separately pinned upstream pure reader, without native imports."""
    path = Path(path)
    if path.stat().st_size > maximum:
        raise ValueError("profile exceeds recorder limit")
    spec = importlib.util.spec_from_file_location("pinned_vmprof", reader_path)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    with path.open("rb") as stream:
        state = reader._read_prof(stream)
    result = sample_summary(
        state.profiles, dict(state.virtual_ips), reader.AssemblerCode
    )
    result.update(
        file=str(path),
        samples=state.n_samples,
        period_us=state.period,
        format=state.version,
    )
    return result


def sample_summary(profiles, symbols, assembler_type):
    """Count logical PyPy frames, retaining compiled-code coverage separately.

    VMProf emits an AssemblerCode PC after a JittedCode logical frame. The PC
    is not another Python function; counting it as the leaf hides the logical
    caller. Inlined operations inside a JIT scope remain unattributed.
    """
    leaf, inclusive = collections.Counter(), collections.Counter()
    weighted, compiled = 0.0, 0.0
    for trace, weight, _, _ in profiles:
        names = []
        has_assembler = False
        for pc in trace:
            if isinstance(pc, assembler_type):
                has_assembler = True
                continue
            names.append(symbols.get(int(pc), "unknown:" + type(pc).__name__))
        if not names:
            continue
        weighted += weight
        if has_assembler:
            compiled += weight
        leaf[names[-1]] += weight
        for name in set(names):
            inclusive[name] += weight
    if not weighted:
        return dict(weighted_samples=0, leaf=[], inclusive=[])
    unknown = sum(v for k, v in leaf.items() if k.startswith("unknown:"))
    return dict(
        weighted_samples=weighted,
        jit_compiled_percent=100 * compiled / weighted,
        unmapped_logical_leaf_percent=100 * unknown / weighted,
        leaf=[
            dict(function=name, percent=100 * count / weighted)
            for name, count in leaf.most_common(25)
        ],
        inclusive=[
            dict(function=name, percent=100 * count / weighted)
            for name, count in inclusive.most_common(25)
        ],
    )
