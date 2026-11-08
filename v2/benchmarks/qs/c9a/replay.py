"""Checked finite replay of complete eligible-atom and boundary journals."""

import time
from collections import defaultdict
from dataclasses import asdict, replace

from v2.benchmarks.infrastructure.resource_control import (
    RecordingBudget as Budget,
)

from .... import qs
from ....qs.polynomial import polynomial_roots
from .configuration import ARMS, configuration, decode
from .journal import records


def replay(
    path,
    arm,
    *,
    seconds=180,
    work_limit=10**13,
    require_mature=True,
    combined=False,
    extended=False,
):
    """Replay G/D/S/P on a fixed admission population, not an I/E experiment.

    Validate the entire stream before publication. Read it again while
    replaying and require the same terminal hash, including after an early
    factor. No discarded source atoms or source ownership are inferred.
    """
    if arm not in ARMS:
        raise ValueError("unknown replay policy")
    started, cpu = time.monotonic(), time.process_time()
    budget = Budget(
        work_limit=work_limit, seconds=seconds, cpu_seconds=seconds
    )
    header, terminal = None, None
    for record in records(path, budget=budget):
        if record["kind"] == "header":
            header = record["data"]
        terminal = record
    if header is None or header.get("format") not in (
        "c9a-atoms-v1",
        "c9a-atoms-v2",
    ):
        raise ValueError("unsupported journal format")
    if (
        require_mature
        and terminal["data"]["result"].get("reason") != "factor_found"
    ):
        raise ValueError(
            "a partial observation is not a mature complete-run stream"
        )
    cell = header["cell"]
    source = decode(cell["configuration"], qs)
    if extended and arm == "C0":
        if getattr(source.collector, "extension_policy", None):
            raise ValueError(
                "a DLP incumbent cannot replay "
                "an expanded admission population"
            )
        extended = combined = False
    if extended:
        from .extended import configuration as expanded_configuration

        if arm not in ("C0", "C1"):
            raise ValueError(
                "expanded grouped diagnosis requires "
                "a separately identified literal policy"
            )
        candidate = (
            source
            if getattr(source.collector, "extension_policy", None)
            else expanded_configuration(
                cell["digits"],
                "C1",
                admission_base_bound=cell["resolved"]["bound"],
                incumbent=source,
            )
        )
    elif combined:
        from .combined import configuration as combined_configuration

        candidate = combined_configuration(
            cell["digits"],
            "C0" if arm == "C0" else "C1",
            admission_base_bound=cell["resolved"]["bound"],
        )
        changes = {
            "noP": dict(preparation_policy="fixed"),
            "noG": dict(cycle_policy="forest"),
            "noD": dict(duplicate_policy="identity"),
            "noS": dict(record_policy="objects"),
        }
        candidate = replace(
            candidate,
            collector=replace(candidate.collector, **changes.get(arm, {})),
        )
    else:
        candidate = configuration(
            cell["digits"], arm, admission_base_bound=cell["resolved"]["bound"]
        )
    # Admission is already fixed by the source stream. Preserve its exact
    # original bounds for independent verification; change only replayable
    # mechanisms. This avoids presenting missing I/E atoms as an ablation.
    cls = (
        qs.BoundedYieldSieveConfig
        if extended
        and getattr(source.collector, "extension_policy", None)
        == "extended-v2"
        else qs.ExtendedYieldSieveConfig
        if extended
        else qs.CombinedYieldSieveConfig
        if combined and arm != "C0"
        else qs.YieldSieveConfig
    )
    # Encoded admission bounds belong to the source stream. Only the chosen
    # structural policies change; recovery cannot invent absent input atoms.
    values = {
        key: value
        for key, value in asdict(source.collector).items()
        if key in cls.__dataclass_fields__
    }
    values.update(
        {
            key: getattr(candidate.collector, key, default)
            for key, default in (
                ("preparation_policy", "fixed"),
                ("cycle_policy", "forest"),
                ("cycle_search_edges", 4096),
                ("duplicate_policy", "identity"),
                ("duplicate_bytes", 16 * 2**20),
                ("record_policy", "objects"),
            )
        }
    )
    values["yield_policy"] = (
        "yield-v2"
        if issubclass(cls, qs.CombinedYieldSieveConfig)
        else "yield-v1"
    )
    if issubclass(cls, qs.CombinedYieldSieveConfig):
        for key in (
            "recovery_policy",
            "split_policy",
            "batch_policy",
            "batch_size",
            "checked_bytes",
            "filter_policy",
            "extraction_policy",
        ):
            values[key] = getattr(candidate.collector, key)
    if extended:
        from ....qs.extended_policy import BOUNDED_FIELDS, EXTENSION_FIELDS

        values.update(
            {
                key: getattr(candidate.collector, key)
                for key in EXTENSION_FIELDS
            }
        )
        if cls is qs.BoundedYieldSieveConfig:
            values.update(
                {
                    key: getattr(candidate.collector, key)
                    for key in BOUNDED_FIELDS
                }
            )
    collector = cls(**values)
    built = qs.build_factor_base(
        cell["n"],
        bound=cell["resolved"]["bound"],
        multiplier=cell["resolved"]["multiplier"],
        budget=budget,
    )
    if built.divisor is not None or built.factor_base is None:
        raise ValueError("replay source lacks its original factor base")
    base, job, divisor = built.factor_base, None, None
    counts, stopped, last = {}, False, None
    for record in records(path, budget=budget):
        kind, data = record["kind"], record["data"]
        last = record
        if kind == "base":
            if data != dict(
                n=int(base.n),
                multiplier=base.multiplier,
                bound=base.bound,
                primes=list(base.primes),
            ):
                raise ValueError("replay factor base changed")
        if stopped:
            continue
        if kind == "window":
            polynomial = qs.Polynomial(
                base.n, base.multiplier, *data["polynomial"]
            )
            if job is None:
                job = qs.QSJob(
                    polynomial,
                    base,
                    data["lo"],
                    data["hi"],
                    config=replace(collector),
                    budget=budget,
                    weight_two=source.weight_two,
                    row_excess=source.row_excess,
                    batch_width=source.batch_width,
                    filter_row_growth=source.filter_row_growth,
                    tested_dependencies=source.tested_dependencies,
                    collection_continues=source.paced_collection,
                )
            else:
                roots = tuple(
                    polynomial_roots(polynomial, base, entry, budget=budget)
                    for entry in base.entries
                )
                job.collector.set_polynomial(
                    polynomial, roots, retain_relations=source.shared_relations
                )
                job.lo = job.next_position = data["lo"]
                job.hi = data["hi"]
                job.final_solve_done = False
                if not source.shared_relations:
                    job.last_count = job.last_solved_count = -1
                    job.zero_checked_count = 0
                    if job._pair_index is not None:
                        job._pair_index.clear()
                    if job._useful_target is not None:
                        job._useful_target = type(job._useful_target)()
                    job._full_view = job._retry_full_view = (
                        job._selected_subset
                    ) = False
        elif kind == "atom_begin":
            if job is None:
                raise ValueError("atom precedes polynomial activation")
            polynomial = qs.Polynomial(
                base.n, base.multiplier, *data["polynomial"]
            )
            if polynomial != job.collector.polynomial:
                raise ValueError("atom is outside its recorded polynomial")
            if data.get("residual_contract") == "extended-v1":
                from ....qs.residual_relations import ResidualAtom

                atom = ResidualAtom(
                    polynomial,
                    data["position"],
                    data["sign"],
                    tuple(tuple(pair) for pair in data["exponents"]),
                    primes=tuple(data["large_primes"]),
                    opaque=data["opaque"],
                )
            else:
                atom = qs.AtomicRelation(
                    polynomial,
                    data["position"],
                    data["sign"],
                    tuple(tuple(pair) for pair in data["exponents"]),
                    data["residual"],
                    large_primes=tuple(data["large_primes"]),
                )
            if atom.relation_id != data["id"]:
                raise ValueError("original replay atom identity changed")
            qs.verify_atomic(
                atom,
                base,
                residual_bound=collector.residual_bound,
                large_prime_bound=collector.large_prime_bound,
                large_product_bound=collector.large_product_bound,
                small_prime_bound=collector.small_prime_bound,
                budget=budget,
                residual_arity=getattr(collector, "residual_arity", 2),
                opaque_pairs=getattr(collector, "opaque_pairs", False),
            )
            # _admit's counters have a closed schema in ordinary collection.
            counters = defaultdict(int, counts)
            refusal = job.collector._admit(atom, counters)
            counts = dict(counters)
            if refusal is not None:
                job.storage_reason = refusal
                divisor = job._solve(True)
                stopped = True
        elif kind == "solve_begin":
            if job is None or data["hi"] != job.hi:
                raise ValueError("solve boundary outside its window")
            job.next_position = data["next_position"]
            if not job.lo <= job.next_position <= job.hi:
                raise ValueError("invalid replay traversal boundary")
            divisor = job._solve(data["final"])
            stopped = divisor is not None
    if last != terminal:
        raise ValueError("journal changed between validation and replay")
    if divisor is not None and not (
        1 < divisor < base.n and base.n % divisor == 0
    ):
        raise ArithmeticError("invalid replay divisor")
    return dict(
        arm=arm,
        source_cell=cell,
        divisor=divisor,
        seconds=time.monotonic() - started,
        cpu_seconds=time.process_time() - cpu,
        work=budget.used,
        stats={} if job is None else job.stats,
        collector_counts=counts,
        stream_terminal=terminal["previous"],
        interpretation="fixed eligible population; no additional I/E yield",
        storage_reason=None if job is None else job.storage_reason,
    )
