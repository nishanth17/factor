"""Compare pinned libraries through the same external C9a timing harness."""

import argparse
import hashlib
import importlib
import importlib.util
import json
import math
import os
import resource
import statistics
import subprocess
import sys
import time
import zipfile
from pathlib import Path

from v2.benchmarks.infrastructure.resource_control import (
    recording_budget,
)

from . import worker
from .configuration import decode
from .journal import canonical
from .launcher import collect_output, host_admission, validate_result


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def archive_library(root, commit, directory):
    """Archive committed runtime files without overlaying another library."""
    root, directory = Path(root), Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    names = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", commit, "v2"],
        cwd=root,
        text=True,
    ).splitlines()
    names = [
        name
        for name in names
        if name.endswith(".py")
        and not name.startswith(("v2/benchmarks/", "v2/tests/"))
    ]
    archive = directory / "library.zip"
    if archive.exists():
        raise ValueError("source archive identity already exists")
    subprocess.run(
        [
            "git",
            "archive",
            "--format=zip",
            "--output=" + str(archive),
            commit,
            *names,
        ],
        cwd=root,
        check=True,
    )
    files = {}
    with zipfile.ZipFile(archive) as zipped:
        for item in zipped.infolist():
            path = Path(item.filename)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("unsafe source archive path")
            if item.is_dir():
                continue
            data = zipped.read(item)
            files[item.filename] = hashlib.sha256(data).hexdigest()
            target = directory / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
    return dict(
        commit=commit,
        root=str(directory.resolve()),
        files=files,
        archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
    )


def load_library(snapshot):
    """Validate runtime bytes before importing a private source package."""
    if snapshot is None:
        from .... import qs
        from ....execution import budget

        return qs, budget
    root = Path(snapshot["root"])
    for name, expected in snapshot["files"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise ValueError("pinned runtime source changed")
    name = "c9a_pinned_" + snapshot["commit"]
    if name not in sys.modules:
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


def prefix_cost(qs, budgets, cell, *, blocks=128, seconds=20):
    """Include setup, collection, verification and store publication cost."""
    budget = recording_budget(
        budgets.Budget, work_limit=10**14, seconds=seconds, cpu_seconds=seconds
    )
    start, cpu = time.monotonic(), time.process_time()
    job = qs.SIQSJob(
        cell["n"],
        seed=cell["seed"],
        config=decode(cell["configuration"], qs),
        budget=budget,
    )
    result = job.run(max_blocks=blocks)
    wall, spent = time.monotonic() - start, time.process_time() - cpu
    if result.reason not in ("paused", "factor_found"):
        raise RuntimeError("fixed-work prefix failed: " + result.reason)
    worker.checked_split(cell["n"], result.divisor)
    collector = job.engine.collector
    atoms = [
        (
            atom.relation_id,
            atom.sign,
            tuple(atom.exponents),
            atom.residual_primes,
        )
        for atom in collector._atoms.values()
    ]
    rows = [
        (
            getattr(row, "atom_ids", (getattr(row, "relation_id", None),)),
            row.sign,
            tuple(row.exponents),
        )
        for row in collector._rows
    ]
    signature = digest(
        dict(
            atoms=atoms,
            rows=rows,
            position=result.next_position,
            divisor=result.divisor,
        )
    )
    transfers = getattr(collector, "_power_transfers", None)
    hit_words = getattr(collector, "_word_hits", None)
    return dict(
        wall=wall,
        cpu=spent,
        work=budget.used,
        signature=signature,
        atoms=len(atoms),
        rows=len(rows),
        positions=job.engine.stats["collected_positions"],
        workspace=collector._workspace,
        helper_bytes=getattr(collector, "_reuse_workspace", 0),
        transfer_builds=0 if transfers is None else transfers.builds,
        transfer_replays=0 if transfers is None else transfers.replays,
        transfer_work=0 if transfers is None else transfers.extra_work,
        transfer_reserved=0 if transfers is None else transfers.used,
        hit_words_work=0 if hit_words is None else hit_words.work,
    )


def measure_rotating(
    functions, groups, *, deadline, samples=9, warm_seconds=3
):
    """Validate every result; retain a fixed rotating order without retries."""
    goldens, warm, observations = {}, {}, []

    def run(label):
        if time.monotonic() >= deadline:
            raise RuntimeError("component measurement deadline")
        row = functions[label]()
        group = groups[label]
        expected = goldens.setdefault(group, row["signature"])
        if row["signature"] != expected:
            raise ArithmeticError("matched component output changed: " + label)
        return dict(mode=label, **row)

    for label in functions:
        wall = cpu = 0
        count = 0
        while min(wall, cpu) < warm_seconds:
            row = run(label)
            wall += row["wall"]
            cpu += row["cpu"]
            count += 1
        warm[label] = dict(wall=wall, cpu=cpu, iterations=count)
    labels = list(functions)
    for index in range(samples):
        order = labels[index % len(labels) :] + labels[: index % len(labels)]
        if index % 2:
            order.reverse()
        for label in order:
            observations.append(dict(sample=index, **run(label)))
    summaries = {}
    for label in labels:
        values = [row["wall"] for row in observations if row["mode"] == label]
        median = statistics.median(values)
        summaries[label] = dict(
            median=median,
            spread=(max(values) - min(values)) / median,
            samples=values,
            stable=(max(values) - min(values)) <= 0.15 * median,
        )
    return dict(warm=warm, observations=observations, summaries=summaries)


def proof_records(qs, base, atoms):
    """Translate verified real records into each library's recovery types."""
    indexed = importlib.import_module(qs.__name__ + ".indexed_relations")
    residual = importlib.import_module(qs.__name__ + ".residual_relations")
    polynomials, records = {}, []
    for source in atoms:
        polynomial = source.polynomial
        fields = (
            polynomial.n,
            polynomial.multiplier,
            polynomial.a,
            polynomial.b,
            polynomial.square_coefficient,
        )
        if fields not in polynomials:
            polynomials[fields] = qs.Polynomial(*fields)
        factors = indexed.IndexedFactors(
            tuple((base._columns[p], e) for p, e in source.exponents),
            base.primes,
        )
        arguments = polynomials[fields], source.position, source.sign, factors
        opaque = getattr(source, "opaque", None)
        record = (
            residual.ResidualAtom(
                *arguments, primes=source.large_primes, opaque=opaque
            )
            if source.large_primes or opaque is not None
            else qs.AtomicRelation(*arguments)
        )
        records.append(record)
    return tuple(records)


def proof_cost(qs, budgets, base, config, records):
    """Time proof, packing, publication and repeated checked reads."""
    indexed = importlib.import_module(qs.__name__ + ".indexed_relations")
    support = importlib.import_module(qs.__name__ + ".yield_support")
    polling = importlib.import_module(
        budgets.__name__.rsplit(".", 1)[0] + ".work_budget"
    )
    budget = polling.PollingBudget(
        recording_budget(
            budgets.Budget,
            work_limit=10**14,
            seconds=20,
            cpu_seconds=20,
        )
    )
    start, cpu = time.monotonic(), time.process_time()
    checked = support.CheckedStore(base, config)
    arena = indexed.RecordArena(base, config=config, budget=budget)
    store = {}
    for record in records:
        checked.atom(record, budget)
        key, divisor = checked.congruence(record, budget)
        if divisor is not None or key is None:
            raise ArithmeticError(
                "real proof replay found a nonunit correction"
            )
        packed, _, slot = arena.plan(record, None, budget)
        arena.commit(packed, slot)
        store[packed.relation_id] = packed
        checked.remember(packed, key=key)
        if hasattr(checked, "retire_proof"):
            checked.retire_proof(record)
    for _ in range(2):
        for packed in store.values():
            checked.atom(packed, budget)
            key, divisor = checked.congruence(packed, budget)
            if divisor is not None or not checked.duplicate(
                key, store, budget
            ):
                raise ArithmeticError(
                    "published canonical proof failed its recheck"
                )
    budget.consume(0)
    wall, spent = time.monotonic() - start, time.process_time() - cpu
    actual = [(a.relation_id, tuple(a.exponents)) for a in store.values()]
    expected = [(a.relation_id, tuple(a.exponents)) for a in records]
    if actual != expected:
        raise ArithmeticError("packed proof payload changed")
    raw_ids = {id(a) for a in records}
    return dict(
        wall=wall,
        cpu=spent,
        signature=digest(actual),
        work=budget.used,
        records=len(records),
        slab_bytes=arena.allocated_bytes,
        proof_records=len(checked.records),
        raw_proofs=sum(identity in raw_ids for identity in checked.records),
        proof_reserved=checked.used,
        canonical_reserved=checked.canonical_used,
    )


def launch_pinned(cell, fixture, snapshot, directory, name):
    """Share the worker, RSS observer, lifecycle and validation rules."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    host = host_admission(cell["rss_bytes"])
    request = directory / (name + ".request.json")
    with request.open("xb") as stream:
        stream.write(canonical(dict(cell=cell, snapshot=snapshot)))
    allowance = cell["seconds"] + cell["warm_seconds"] + cell["margin"]

    def limits():
        resource.setrlimit(
            resource.RLIMIT_CPU, (math.ceil(allowance), math.ceil(allowance))
        )
        resource.setrlimit(resource.RLIMIT_FSIZE, (257 * 2**20, 257 * 2**20))

    started = time.monotonic()
    with (
        (directory / (name + ".stdout.json")).open("xb") as out,
        (directory / (name + ".stderr.txt")).open("xb") as err,
    ):
        process = subprocess.Popen(
            [sys.executable, "-B", "-m", __name__, "--worker", str(request)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
            preexec_fn=limits,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        )
        reason = collect_output(
            process, out, err, started + allowance, cell["rss_bytes"]
        )
    if reason is not None or process.returncode:
        # Preserve all outputs and the reserved charge; no automatic retry.
        raise RuntimeError(
            "pinned worker failed: " + str((reason, process.returncode))
        )
    result = json.loads((directory / (name + ".stdout.json")).read_text())
    intended = (
        Path(snapshot["root"]) / "v2"
        if snapshot
        else Path(__file__).resolve().parents[4] / "v2"
    )
    if result["library_root"] != str(intended.resolve()):
        raise ValueError("worker loaded the wrong runtime")
    result = validate_result(result, fixture, cell)
    result.update(
        host=host,
        source_commit=None if snapshot is None else snapshot["commit"],
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    args = parser.parse_args()
    request = json.loads(args.worker.read_text())
    worker.candidate, worker.candidate_budget = load_library(
        request["snapshot"]
    )
    result = worker.run(request["cell"])
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
