"""Compact checked SIQS stores and replayable progress; no matrix dumps."""

import hashlib
import json
from dataclasses import asdict
from functools import partial

from ..common import arithmetic, utils
from ..execution.resources import decode_checkpoint, encode_checkpoint
from .assignment_stream import assignment_identity
from .extended_policy import (
    BOUNDED_FIELDS,
    EXTENSION_FIELDS,
    BoundedYieldSieveConfig,
    ExtendedYieldSieveConfig,
)
from .extraction import DependencyExtractor, _prepare_relations
from .families import (
    PolynomialFamily,
    _checked_resources,
    _checksum,
    _identity,
)
from .linear_algebra import DependencySolver, filter_matrix
from .pipeline import QSJob
from .polynomial import Polynomial, polynomial_roots
from .relations import (
    AtomicRelation,
    CombinedRelation,
    _combination_workspace,
    combined_storage_reserve,
    verify_atomic,
    verify_combined,
)
from .sieve_collector import (
    CombinedYieldSieveConfig,
    DoubleLargeSieveConfig,
    GrowingDoubleLargeSieveConfig,
    SieveCollector,
    SieveConfig,
    YieldSieveConfig,
)
from .storage import require_published

VERSION = 3
MAX_BLOB_BYTES = 64 * 1024 * 1024


def _solver_digest(solver, *, encoding="decimal-v1"):
    """Retain legacy fingerprints; new masks use bounded hexadecimal text."""
    state = [
        solver.next_row,
        solver.xors,
        solver.pending,
        sorted(solver.pivots.items()),
        solver.dependencies,
    ]
    if encoding == "hex-v1":
        digest = hashlib.sha256()

        def feed(value):
            if arithmetic.is_integer(value):
                digest.update(('"' + hex(value) + '"').encode())
            elif value is None:
                digest.update(b"null")
            else:
                digest.update(b"[")
                for index, item in enumerate(value):
                    if index:
                        digest.update(b",")
                    feed(item)
                digest.update(b"]")

        feed(state)
        return digest.hexdigest()

    if encoding != "decimal-v1":
        raise ValueError("unknown SIQS solver digest encoding")
    return _checksum(state)


def _store(collector):
    require_published(collector)
    if collector._graph is not None:
        from .large_prime_checkpoint import pack_store

        return pack_store(collector)
    atoms = list(collector._atoms.values())
    indices = {atom.relation_id: i for i, atom in enumerate(atoms)}
    polynomials, lookup, encoded = [], {}, []

    for atom in atoms:
        key = (
            atom.polynomial.a,
            atom.polynomial.b,
            atom.polynomial.square_coefficient,
        )
        if key not in lookup:
            lookup[key] = len(polynomials)
            polynomials.append(key)
        encoded.append(
            [
                lookup[key],
                atom.position,
                atom.sign,
                atom.exponents,
                atom.residual,
            ]
        )

    combined = [
        [
            tuple(indices[i] for i in item.atom_ids),
            item.u,
            item.sign,
            item.exponents,
            item.square_correction,
        ]
        for item in collector._combined
    ]
    rows = collector._full + collector._combined
    # Solver masks use mixed admission order, while storage encodes full
    # and combined payloads separately. Retain the mapping between them.
    row_indices = {id(row): index for index, row in enumerate(rows)}
    return dict(
        row_order=[row_indices[id(row)] for row in collector._rows],
        polynomials=polynomials,
        atoms=encoded,
        full=[indices[item.relation_id] for item in collector._full],
        combined=combined,
        pending=[[r, indices[i]] for r, i in collector._pending.items()],
    )


def _arithmetic_version(values):
    """Identify coupled arithmetic independently of the collector schema."""
    if (
        values.get("root_update_policy", "standard") != "standard"
        or values.get("assignment_policy") == "windowed_pair"
        or values.get("multiplier_policy") == "parity2_norm"
    ):
        return 10
    marking = values["collector"].get("marking", "")
    if (
        values.get("multiplier_policy") == "parity2_wide"
        or values.get("assignment_policy") == "paired_flyer"
        or (isinstance(marking, str) and "parity_two" in marking.split("+"))
    ):
        return 9
    if (
        values.get("polynomial_policy", "standard") != "standard"
        or values.get("multiplier_policy", "standard") != "standard"
    ):
        return 8
    return 0


@encode_checkpoint
def pack_job(job):
    """Snapshot immutable provenance once, with bounded encoded output."""
    config_values = asdict(job.config)
    arithmetic_version = _arithmetic_version(config_values)
    version = (
        14
        if isinstance(job.config.collector, BoundedYieldSieveConfig)
        else 13
        if isinstance(job.config.collector, ExtendedYieldSieveConfig)
        else 12
        if isinstance(job.config.collector, CombinedYieldSieveConfig)
        else 11
        if isinstance(job.config.collector, YieldSieveConfig)
        else arithmetic_version
        if arithmetic_version
        else 7
        if job.config.sign_policy != "canonical"
        or job.config.collector.sparse_queue_bytes
        else 6
        if job.config.solve_policy != "window"
        else 5
        if isinstance(job.config.collector, GrowingDoubleLargeSieveConfig)
        else 4
        if isinstance(job.config.collector, DoubleLargeSieveConfig)
        else VERSION
    )
    engine = job.engine
    store = _store(engine.collector) if engine is not None else None
    progress = None
    if engine is not None:
        solver = engine.solver
        progress = dict(
            polynomial=[
                engine.collector.polynomial.a,
                engine.collector.polynomial.b,
                engine.collector.polynomial.square_coefficient,
            ],
            lo=engine.lo,
            hi=engine.hi,
            next_position=engine.next_position,
            final_solve_done=engine.final_solve_done,
            storage_reason=engine.storage_reason,
            storage_solve_done=engine.storage_solve_done,
            last_count=engine.last_count,
            last_solved_count=engine.last_solved_count,
            stats=engine.stats,
            workspace_peak=engine._workspace_peak,
            solver=None
            if solver is None
            else dict(
                actions=solver.next_row + solver.xors,
                digest=_solver_digest(solver, encoding="hex-v1"),
                digest_encoding="hex-v1",
                pending=solver.pending is not None,
                extract_index=None
                if engine.extractor is None
                else engine.extractor.next_dependency,
            ),
        )
        if version >= 6:
            progress["zero_checked_count"] = engine.zero_checked_count
        if version in (13, 14) and solver is not None:
            progress["solver"]["incremental_modes"] = [
                solver.matrix.stats.get("incremental_survivors", False),
                solver.matrix.stats.get("incremental_core", False),
            ]
        if version in (11, 12, 13, 14):
            progress["yield_progress"] = {
                "pair_checked_count": 0
                if engine._pair_index is None
                else engine._pair_index.checked_count
            }
            if version >= 12:
                progress["yield_progress"].update(
                    active_target=None
                    if engine._useful_target is None
                    else engine._useful_target.state(),
                    solve_view=[
                        engine._full_view,
                        engine._retry_full_view,
                        engine._selected_subset,
                    ],
                )
                if progress["solver"] is not None:
                    progress["solver"]["light_rows"] = engine._selected_subset

    payload = dict(
        version=version,
        backend=arithmetic.get_backend(job.config.backend).identity,
        n=job.n,
        seed=job.seed,
        config=config_values,
        multiplier=job.multiplier,
        base_identity=None if job.base is None else _identity(job.base),
        assignment_identity=None
        if job.assignments is None
        else assignment_identity(job.assignments),
        family_index=job.family_index,
        coefficient_cursor=job.coefficient_cursor,
        gray_index=None if job.family is None else job.family.next_index,
        pending=None
        if job.pending_step is None
        else [
            job.pending_step[0].a,
            job.pending_step[0].b,
            job.pending_step[0].square_coefficient,
        ],
        epoch=job.epoch,
        half_width=job.half_width,
        stalled=job.stalled,
        active=job.active,
        seen=sorted(job.seen),
        finished_reason=job.finished_reason,
        divisor=job.divisor,
        poly_start_rows=job.poly_start_rows,
        stats=job.stats,
        engine=progress,
        store=store,
        store_identity=_checksum(store),
    )
    if arithmetic_version == 10:
        payload["coupled_semantics"] = "c10a01-coupled-v4"
    else:
        payload["config"].pop("root_update_policy")
    if arithmetic_version == 9:
        payload["coupled_semantics"] = "c10a01-coupled-v3"
    if version >= 6:
        payload["pending_final_reason"] = job.pending_final_reason
    else:
        # Keep legacy checkpoint configuration and schema identities intact.
        payload["config"].pop("solve_policy")
    if arithmetic_version < 8:
        payload["config"].pop("polynomial_policy")
        payload["config"].pop("multiplier_policy")
    if version < 7:
        payload["config"].pop("sign_policy")
        payload["config"]["collector"].pop("sparse_queue_bytes")
    encoder = json.JSONEncoder(
        default=arithmetic.json_integer, sort_keys=True, separators=(",", ":")
    )
    parts, size = [], 0
    for part in encoder.iterencode(payload):
        size += len(part.encode())
        parts.append(part)

    blob = "".join(parts)
    digest = hashlib.sha256(blob.encode()).hexdigest()
    resources = _checked_resources(
        dict(
            work_used=job.budget.used,
            wall_used=job.budget.wall_used,
            cpu_used=job.budget.cpu_used,
        )
    )
    return dict(
        version=version,
        blob=blob,
        sha256=digest,
        resources=resources,
        resources_sha256=_checksum(resources),
    )


def _restore_store(payload, collector, budget):
    """Reserve decoded provenance and verify every atom/combination afresh."""
    if collector._graph is not None:
        from .large_prime_checkpoint import restore_store

        return restore_store(payload, collector, budget)
    if not isinstance(payload, dict) or set(payload) - {"row_order"} != {
        "polynomials",
        "atoms",
        "full",
        "combined",
        "pending",
    }:
        raise ValueError("invalid SIQS relation store shape")

    config, base = collector.config, collector.factor_base

    for name in ("polynomials", "atoms", "full", "combined", "pending"):
        if not isinstance(payload[name], list):
            raise ValueError("invalid SIQS store sequence")
    if (
        config.max_partials is not None
        and len(payload["pending"]) > config.max_partials
    ):
        raise ValueError("saved partial FIFO exceeds selected capacity")

    row_count = len(payload["full"]) + len(payload["combined"])
    order = (
        payload["row_order"]
        if "row_order" in payload
        else list(range(row_count))
    )
    if not isinstance(order, list) or len(order) != row_count:
        raise ValueError("invalid checkpoint mixed row order")
    for index in order:
        utils.require_integer(index, "matrix row index", 0)
        if index >= row_count:
            raise ValueError("checkpoint row order refers outside store")
    if len(set(order)) != row_count:
        raise ValueError("checkpoint row order repeats a row")
    reserves = []

    for record in payload["atoms"]:
        if not isinstance(record, list) or len(record) != 5:
            raise ValueError("invalid SIQS checkpoint atom")
        _, position, _, exponents, _ = record
        utils.require_integer(position, "position")
        if (
            abs(position).bit_length() > 4096
            or not isinstance(exponents, list)
            or len(exponents) > len(base.entries)
        ):
            raise ValueError("checkpoint atom size exceeds its bound")

        reserves.append(
            4096
            + 256 * len(exponents)
            + 16 * (abs(position).bit_length() + base.n_prime.bit_length())
        )

    combined_reserve = 0

    for record in payload["combined"]:
        if not isinstance(record, list) or len(record) != 5:
            raise ValueError("invalid checkpoint match")
        indices = record[0]
        if not isinstance(indices, list) or len(indices) != 2:
            raise ValueError("checkpoint match needs two atoms")
        for index in indices:
            utils.require_integer(index, "atom index", 0)
            if index >= len(payload["atoms"]):
                raise ValueError("missing checkpoint provenance")
        combined_reserve += combined_storage_reserve(
            sum(len(payload["atoms"][i][3]) for i in indices),
            base.n.bit_length(),
        )

    (collector._workspace + sum(reserves) + combined_reserve)
    polynomials = [
        Polynomial(base.n, base.multiplier, *values)
        for values in payload["polynomials"]
    ]
    if len(set(p.identity for p in polynomials)) != len(polynomials):
        raise ValueError("duplicate checkpoint polynomial")
    atoms = []

    for index, (poly_index, position, sign, exponents, residual) in enumerate(
        payload["atoms"]
    ):
        utils.require_integer(poly_index, "polynomial index", 0)
        if poly_index >= len(polynomials):
            raise ValueError("checkpoint atom polynomial is missing")
        atom = AtomicRelation(
            polynomials[poly_index],
            position,
            sign,
            tuple(tuple(pair) for pair in exponents),
            residual,
        )
        verify_atomic(
            atom, base, residual_bound=config.residual_bound, budget=budget
        )
        if atom.relation_id in collector._atoms:
            raise ValueError("duplicate checkpoint atomic identity")
        collector._atoms[atom.relation_id] = atom
        collector._atom_bytes[atom.relation_id] = reserves[index]
        atoms.append(atom)

    used = set()

    def take(index):
        utils.require_integer(index, "atom index", 0)
        if index >= len(atoms) or index in used:
            raise ValueError("missing or reused checkpoint provenance")
        used.add(index)
        return atoms[index]

    for index in payload["full"]:
        atom = take(index)
        if atom.residual != 1:
            raise ValueError("partial checkpoint atom is marked full")
        collector._full.append(atom)

    for indices, u, sign, exponents, correction in payload["combined"]:
        if len(indices) != 2:
            raise ValueError("checkpoint match needs two atoms")
        selected = [take(i) for i in indices]
        item = CombinedRelation(
            tuple(a.relation_id for a in selected),
            arithmetic.backend_for(base.n).integer(u),
            sign,
            tuple(tuple(pair) for pair in exponents),
            arithmetic.backend_for(base.n).integer(correction),
        )
        scratch = _combination_workspace(selected, base) - base.workspace_bytes
        collector._scratch_peak_bytes = max(
            collector._scratch_peak_bytes, scratch
        )
        verify_combined(
            item,
            base,
            collector._atoms,
            budget=budget,
        )
        collector._combined.append(item)
        last = atoms[max(indices)].relation_id
        collector._atom_bytes[last] += combined_storage_reserve(
            sum(len(a.exponents) for a in selected), base.n.bit_length()
        )

    for residual, index in payload["pending"]:
        atom = take(index)
        if (
            residual != atom.residual
            or residual == 1
            or residual in collector._pending
        ):
            raise ValueError("invalid checkpoint partial matching state")

        collector._pending[residual] = atom.relation_id

    if len(used) != len(atoms):
        raise ValueError("checkpoint has unreferenced atomic provenance")
    rows = collector._full + collector._combined
    collector._rows.extend(rows[index] for index in order)
    collector._workspace += sum(collector._atom_bytes.values())


@decode_checkpoint
def restore_job(checkpoint, *, budget, config=None, allow_extension=False):
    """Verify and rebuild charged roots, stores and solver prefixes."""
    from .siqs import CountedSIQSConfig, SIQSConfig, SIQSJob, make_family

    if not isinstance(checkpoint, dict) or set(checkpoint) != {
        "version",
        "blob",
        "sha256",
        "resources",
        "resources_sha256",
    }:
        raise ValueError("invalid full SIQS checkpoint envelope")

    blob = checkpoint["blob"]
    if type(checkpoint["version"]) is not int or checkpoint["version"] not in (
        1,
        2,
        VERSION,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
    ):
        raise ValueError("unsupported SIQS checkpoint version")

    if hashlib.sha256(blob.encode()).hexdigest() != checkpoint["sha256"]:
        raise ValueError("SIQS checkpoint integrity mismatch")
    resources = _checked_resources(checkpoint["resources"])
    if _checksum(resources) != checkpoint["resources_sha256"]:
        raise ValueError("SIQS checkpoint resource integrity mismatch")
    payload = json.loads(blob)
    if type(payload["version"]) is not int or (
        payload["version"] != checkpoint["version"]
    ):
        raise ValueError("invalid SIQS checkpoint payload version")
    if payload["version"] >= 2 and payload["store"] is not None:
        if not isinstance(payload["store"], dict) or (
            "row_order" not in payload["store"]
        ):
            raise ValueError("mixed-order checkpoint lacks row_order")

    values = dict(payload["config"])
    ordered_fields = {
        "yield_policy",
        "small_prime_bound",
        "preparation_policy",
        "cycle_policy",
        "cycle_search_edges",
        "duplicate_policy",
        "duplicate_bytes",
        "record_policy",
    }
    present = ordered_fields.intersection(values["collector"])
    combined_fields = {
        "recovery_policy",
        "split_policy",
        "batch_policy",
        "batch_size",
        "checked_bytes",
        "filter_policy",
        "extraction_policy",
    }
    if payload["version"] in (11, 12, 13, 14):
        if present != ordered_fields or (
            values["collector"].get("yield_policy")
            != ("yield-v2" if payload["version"] >= 12 else "yield-v1")
            or values["collector"].get("storage_policy") != "grow"
        ):
            raise ValueError("yield checkpoint lacks its policy identity")
    elif present:
        raise ValueError("legacy checkpoint cannot change admission policy")
    if combined_fields.intersection(values["collector"]) != (
        combined_fields if payload["version"] >= 12 else set()
    ):
        raise ValueError("combined checkpoint lacks its policy identity")
    if EXTENSION_FIELDS.intersection(values["collector"]) != (
        EXTENSION_FIELDS if payload["version"] in (13, 14) else set()
    ):
        raise ValueError("expanded checkpoint lacks its policy identity")
    if BOUNDED_FIELDS.intersection(values["collector"]) != (
        BOUNDED_FIELDS if payload["version"] == 14 else set()
    ):
        raise ValueError("bounded checkpoint lacks its policy identity")
    if payload["version"] in (13, 14) and values["collector"].get(
        "extension_policy"
    ) != ("extended-v2" if payload["version"] == 14 else "extended-v1"):
        raise ValueError("bounded checkpoint policy differs from schema")
    if payload["version"] < 7 and (
        (values.get("solve_policy", "window") != "window")
        != (payload["version"] == 6)
    ):
        raise ValueError("SIQS scheduling policy disagrees with schema")
    if payload["version"] < 7 and (
        "sign_policy" in values or "sparse_queue_bytes" in values["collector"]
    ):
        raise ValueError("legacy checkpoint cannot change collection policy")
    if payload["version"] >= 7 and (
        "sign_policy" not in values
        or "sparse_queue_bytes" not in values["collector"]
    ):
        raise ValueError("collection checkpoint lacks its policy identity")
    policies = {"polynomial_policy", "multiplier_policy"}
    # C9a schemas11–13 predate coupled-policy fields. Old checkpoints keep
    # standard arithmetic; new combinations bind the C10a marker separately
    # without losing the collector's indexed/provenance/resume schema.
    arithmetic_version = (
        _arithmetic_version(values)
        if payload["version"] >= 11
        else payload["version"]
    )
    if arithmetic_version < 8 and policies.intersection(values):
        raise ValueError("legacy checkpoint cannot change arithmetic policy")
    if arithmetic_version >= 8 and not policies.issubset(values):
        raise ValueError("parity checkpoint lacks its policy identity")
    introduced_v4 = (
        values.get("root_update_policy", "standard") != "standard"
        or values.get("assignment_policy") == "windowed_pair"
        or values.get("multiplier_policy") == "parity2_norm"
    )
    if arithmetic_version < 10 and (
        "root_update_policy" in values or introduced_v4
    ):
        raise ValueError("legacy checkpoint cannot change v4 coupled policy")
    if arithmetic_version == 10 and (
        "root_update_policy" not in values
        or not introduced_v4
        or payload.get("coupled_semantics") != "c10a01-coupled-v4"
    ):
        raise ValueError("v4 checkpoint lacks its coupled policy identity")
    marking = values["collector"].get("marking", "")
    introduced = (
        values.get("multiplier_policy") == "parity2_wide"
        or values.get("assignment_policy") == "paired_flyer"
        or (isinstance(marking, str) and "parity_two" in marking.split("+"))
    )
    if arithmetic_version < 9 and (
        introduced or "coupled_semantics" in payload
    ):
        raise ValueError("legacy checkpoint cannot change coupled policy")
    if arithmetic_version == 9 and (
        not introduced
        or payload.get("coupled_semantics") != "c10a01-coupled-v3"
    ):
        raise ValueError("coupled checkpoint lacks its policy identity")
    collector_type = (
        BoundedYieldSieveConfig
        if payload["version"] == 14
        else ExtendedYieldSieveConfig
        if payload["version"] == 13
        else CombinedYieldSieveConfig
        if payload["version"] == 12
        else YieldSieveConfig
        if payload["version"] == 11
        else GrowingDoubleLargeSieveConfig
        if payload["version"] == 5
        or payload["version"] >= 6
        and values["collector"].get("storage_policy") == "grow"
        else DoubleLargeSieveConfig
        if payload["version"] == 4
        or payload["version"] >= 6
        and "large_prime_bound" in values["collector"]
        else SieveConfig
    )
    values["collector"] = collector_type(**values["collector"])
    config_type = CountedSIQSConfig if "base_count" in values else SIQSConfig
    saved_config = config_type(**values)
    if type(allow_extension) is not bool:
        raise TypeError("allow_extension must be Boolean")
    extended = config is not None and config != saved_config
    if extended:
        if not allow_extension:
            raise ValueError("SIQS checkpoint configuration mismatch")
        from .capacity import validate_extension

        validate_extension(saved_config, config)
    else:
        config = saved_config

    expected_backend = arithmetic.get_backend(config.backend).identity
    if (
        payload.get(
            "backend", "python-int" if payload["version"] < 3 else None
        )
        != expected_backend
    ):
        raise ValueError("incompatible checkpoint backend")

    budget.restore_charges(
        work=resources["work_used"],
        wall=resources["wall_used"],
        cpu=resources["cpu_used"],
    )

    public_config = config
    dynamic_snapshot = True
    budget.consume(0)
    job = SIQSJob(
        payload["n"], seed=payload["seed"], config=config, budget=budget
    )

    for name, limit in (
        ("epoch", config.growth_steps),
        ("family_index", config.family_count),
        ("stalled", config.polynomial_limit),
        (
            "poly_start_rows",
            len(payload["store"]["row_order"]) if payload["store"] else 0,
        ),
    ):
        value = utils.require_integer(payload[name], name, 0)
        if value > limit:
            raise ValueError("SIQS checkpoint progress exceeds its bound")
        setattr(job, name, value)

    expected_width = min(
        config.max_half_width, config.half_width * (2**job.epoch)
    )
    if payload["half_width"] != expected_width:
        raise ValueError("SIQS checkpoint width disagrees with recovery epoch")
    job.half_width = expected_width
    for name in ("active",):
        if type(payload[name]) is not bool:
            raise ValueError("invalid SIQS checkpoint flag")
        setattr(job, name, payload[name])
    job.coefficient_cursor = utils.require_integer(
        payload.get("coefficient_cursor", 0), "coefficient cursor", 0
    )
    if job.coefficient_cursor.bit_length() > 2048:
        raise ValueError("coefficient cursor exceeds its bit limit")
    if not config.external_coefficients and job.coefficient_cursor:
        raise ValueError("nonexternal job has a coefficient cursor")
    if config.external_coefficients and job.coefficient_cursor:
        latest = payload["pending"]
        if latest is None and payload["engine"] is not None:
            latest = payload["engine"]["polynomial"]
        if (
            latest is None
            or len(latest) != 3
            or latest[2] + 4 != job.coefficient_cursor
            or latest[2] % 4 != 3
            or latest[0] != latest[2] ** 2
        ):
            raise ValueError(
                "external coefficient cursor disagrees with progress"
            )

    if config.streaming and payload["seen"]:
        raise ValueError(
            "streaming checkpoint must not contain a history registry"
        )
    if len(payload["seen"]) > config.polynomial_limit:
        raise ValueError("SIQS checkpoint registry exceeds its bound")
    job.seen = set(tuple(key) for key in payload["seen"])
    if len(job.seen) != len(payload["seen"]):
        raise ValueError("duplicate SIQS checkpoint registry entry")
    for a, b, width in job.seen:
        utils.require_integer(a, "registry A", 1)
        utils.require_integer(b, "registry B", 0)
        utils.require_integer(width, "registry width", 1)
        if (
            a.bit_length() > 4096
            or b > a // 2
            or width > config.max_half_width
        ):
            raise ValueError("invalid SIQS checkpoint polynomial registry")

    job.stats = payload["stats"]
    job.divisor = payload["divisor"]
    if job.divisor is not None and not utils.valid_divisor(job.divisor, job.n):
        raise ValueError("invalid checkpoint divisor")
    allowed = (
        None,
        "families_exhausted",
        "assignment_space_exhausted",
        "coefficient_limit",
        "stalled_yield",
        "trivial_dependency_limit",
        "configuration_refusal",
    )
    if payload["finished_reason"] not in allowed:
        raise ValueError("invalid SIQS checkpoint stop reason")
    job.finished_reason = payload["finished_reason"]
    if payload["version"] >= 6:
        pending = payload["pending_final_reason"]
        if pending not in allowed or (
            pending is not None
            and (
                not config.paced_collection
                or job.finished_reason is not None
                or payload["engine"] is None
            )
        ):
            raise ValueError("invalid pending SIQS collection stop")
        job.pending_final_reason = pending
    job.multiplier = payload["multiplier"]
    if job.multiplier is not None:
        utils.require_integer(job.multiplier, "checkpoint multiplier", 1)
        if (
            config.multiplier and job.multiplier != config.multiplier
        ) or job.multiplier > 1000000:
            raise ValueError(
                "SIQS checkpoint multiplier differs from configuration"
            )

    if payload["base_identity"] is not None:
        job._setup()
        if job.base is None or _identity(job.base) != payload["base_identity"]:
            raise ValueError("SIQS checkpoint base identity mismatch")
        if (
            payload["assignment_identity"] is not None
            and assignment_identity(job.assignments)
            != payload["assignment_identity"]
        ):
            raise ValueError("SIQS checkpoint assignment identity mismatch")

        if payload["assignment_identity"] is None:
            job.assignments = None
    elif payload["engine"] is not None or payload["store"] is not None:
        raise ValueError("SIQS checkpoint store has no factor base")

    gray = payload["gray_index"]
    if gray is not None:
        if job.assignments is None or job.family_index >= len(job.assignments):
            raise ValueError("SIQS checkpoint family index is exhausted")
        family_factory = (
            partial(make_family, config=config)
            if config.sign_policy != "canonical"
            or config.polynomial_policy != "standard"
            or config.root_update_policy != "standard"
            else PolynomialFamily
        )
        job.family = family_factory(
            job.base,
            job.assignments[job.family_index],
            budget=budget,
        )
        utils.require_integer(gray, "Gray index", 0)
        if gray > min(job.family.count, config.gray_limit):
            raise ValueError("SIQS checkpoint Gray index exceeds its family")
        if gray:
            job.family.current = job.family._make_step(gray - 1, None)
        job.family.next_index = gray
        job.stats["root_reconstructions"] += 1

    if payload["pending"] is not None:
        polynomial = Polynomial(job.n, job.multiplier, *payload["pending"])
        if job.family is not None and (
            job.family.current is None
            or polynomial != job.family.current.polynomial
        ):
            raise ValueError(
                "SIQS checkpoint pending polynomial disagrees with family"
            )

        roots = tuple(
            polynomial_roots(polynomial, job.base, e, budget=budget)
            for e in job.base.entries
        )
        job.pending_step = (polynomial, roots)

    saved = payload["engine"]
    if saved is not None:
        polynomial = Polynomial(job.n, job.multiplier, *saved["polynomial"])
        roots = tuple(
            polynomial_roots(polynomial, job.base, e, budget=budget)
            for e in job.base.entries
        )
        extra = config.metadata_reserve + 65536 + 640 * len(job.base.entries)
        extra += 16 * (
            max(2, config.factor_count) * job.base.bound.bit_length()
            + job.n.bit_length()
        )
        from dataclasses import replace

        collector_config = replace(config.collector)
        job.engine = QSJob(
            polynomial,
            job.base,
            saved["lo"],
            saved["hi"],
            config=collector_config,
            budget=budget,
            weight_two=config.weight_two,
            row_excess=config.row_excess,
            batch_width=config.batch_width,
            filter_row_growth=config.filter_row_growth,
            tested_dependencies=config.tested_dependencies,
            collection_continues=config.paced_collection,
            collector_class=partial(SieveCollector, precomputed_roots=roots),
        )
        engine = job.engine
        if _checksum(payload["store"]) != payload["store_identity"]:
            raise ValueError(
                "SIQS checkpoint relation store identity mismatch"
            )
        _restore_store(payload["store"], engine.collector, budget)
        if payload["version"] in (11, 12, 13, 14):
            progress = saved.get("yield_progress")
            expected = {"pair_checked_count"}
            if payload["version"] >= 12:
                expected.update(("active_target", "solve_view"))
            if not isinstance(progress, dict) or set(progress) != expected:
                raise ValueError(
                    "yield checkpoint lacks its progress identity"
                )
            count = progress["pair_checked_count"]
            if engine._pair_index is None:
                if type(count) is not int or count != 0:
                    raise ValueError("disabled duplicate policy has progress")
            else:
                engine._pair_index.rebuild(
                    engine.collector._rows, count, job.base, budget
                )
            if payload["version"] >= 12:
                if engine._useful_target is None:
                    if progress["active_target"] is not None:
                        raise ValueError("disabled useful target has progress")
                else:
                    engine._useful_target.restore(
                        progress["active_target"],
                        engine.collector._rows,
                        engine.collector,
                        budget,
                    )
                view = progress["solve_view"]
                if (
                    not isinstance(view, list)
                    or len(view) != 3
                    or any(type(v) is not bool for v in view)
                ):
                    raise ValueError("invalid solve-view progress")
                (
                    engine._full_view,
                    engine._retry_full_view,
                    engine._selected_subset,
                ) = view
                if engine._full_view != engine._retry_full_view or (
                    engine._selected_subset and engine._full_view
                ):
                    raise ValueError("inconsistent solve-view progress")
        elif "yield_progress" in saved:
            raise ValueError("legacy checkpoint cannot change yield progress")
        if payload["version"] >= 6:
            engine.zero_checked_count = utils.require_integer(
                saved["zero_checked_count"], "zero-row cursor", 0
            )
            if engine.zero_checked_count > len(engine.collector._rows):
                raise ValueError("zero-row cursor exceeds retained rows")
        for name in (
            "next_position",
            "last_count",
            "last_solved_count",
            "workspace_peak",
        ):
            value = utils.require_integer(
                saved[name],
                name,
                -1
                if name.startswith("last")
                else 0
                if name == "workspace_peak"
                else saved["lo"],
            )
            if name == "next_position" and value > saved["hi"]:
                raise ValueError(
                    "SIQS checkpoint block position exceeds window"
                )
            setattr(
                engine,
                "_workspace_peak" if name == "workspace_peak" else name,
                value,
            )

        for name in ("final_solve_done", "storage_solve_done"):
            if type(saved[name]) is not bool:
                raise ValueError("invalid SIQS engine completion flag")
            setattr(engine, name, saved[name])
        engine.storage_reason = saved["storage_reason"]
        if engine.storage_reason not in (None,):
            raise ValueError("invalid SIQS checkpoint storage reason")

        engine.stats = saved["stats"]
        prefix = saved["solver"]
        if prefix is not None:
            incremental_modes = prefix.get("incremental_modes", [False, False])
            if payload["version"] in (13, 14) and (
                not isinstance(incremental_modes, list)
                or len(incremental_modes) != 2
                or any(type(flag) is not bool for flag in incremental_modes)
                or (incremental_modes[1] and not incremental_modes[0])
            ):
                raise ValueError("invalid incremental solver modes")
            if payload["version"] >= 12 and (
                type(prefix.get("light_rows")) is not bool
                or prefix["light_rows"] != engine._selected_subset
            ):
                raise ValueError("solver view differs from yield progress")
            relations = engine.collector.matrix_relations
            prepared = _prepare_relations(
                relations,
                job.base,
                engine.collector._atoms,
                budget=budget,
                retained_workspace_bytes=engine.collector._workspace,
                verification_cache=engine.collector._preparation_cache
                if payload["version"] >= 12
                else None,
            )
            (
                engine.collector._workspace
                + prepared.workspace_bytes
                - prepared.shared_workspace_bytes
            )
            matrix = filter_matrix(
                prepared.rows,
                weight_two=config.weight_two,
                budget=budget,
                merge_policy="fill-v2"
                if getattr(collector_config, "filter_policy", "reference")
                != "reference"
                else "reference",
                light_rows=engine._selected_subset,
                surplus=config.row_excess,
                merge_weight=getattr(collector_config, "merge_weight", 2),
                merge_fill=getattr(collector_config, "merge_fill", 0),
                sparse_incidence=getattr(
                    collector_config, "sparse_incidence", False
                ),
                panel_bits=getattr(collector_config, "solver_panel_bits", 0),
                _incremental=engine.collector._incremental
                if getattr(collector_config, "incremental_bytes", 0)
                and incremental_modes[0]
                else None,
                _component_contractions=incremental_modes[1],
            )
            solver = DependencySolver(matrix, budget=budget)
            actions = utils.require_integer(
                prefix["actions"], "solver actions", 0
            )
            if actions > resources["work_used"] or actions > len(
                matrix.rows
            ) ** 2 + len(matrix.rows):
                raise ValueError(
                    "SIQS checkpoint solver prefix exceeds consumed work"
                )

            # Rebuild from checked rows rather than trusting saved pivots.
            # Replay consumes the resumed allowance and verifies its digest.
            for _ in range(actions):
                if solver.next_row == len(matrix.rows):
                    raise ValueError(
                        "SIQS checkpoint solver prefix exceeds matrix"
                    )
                solver.step()

            if (
                prefix["pending"]
                and solver.pending is None
                and solver.next_row < len(matrix.rows)
            ):
                solver.pending = (
                    matrix.rows[solver.next_row],
                    matrix.masks[solver.next_row],
                )

            if (
                _solver_digest(
                    solver,
                    encoding=prefix.get("digest_encoding", "decimal-v1"),
                )
                != prefix["digest"]
            ):
                raise ValueError("SIQS checkpoint elimination prefix mismatch")

            engine.prepared, engine.solver = prepared, solver
            if prefix["extract_index"] is not None:
                if solver.next_row != len(matrix.rows):
                    raise ValueError(
                        "SIQS checkpoint extraction precedes elimination"
                    )
                index = utils.require_integer(
                    prefix["extract_index"], "extraction index", 0
                )
                if index > len(solver.dependencies):
                    raise ValueError(
                        "SIQS checkpoint extraction index exceeds kernels"
                    )
                engine.extractor = DependencyExtractor(
                    prepared,
                    tuple(solver.dependencies),
                    budget=budget,
                    product_policy=getattr(
                        collector_config, "extraction_policy", "reference"
                    ),
                )
                # Completed trials are independently replayed and charged.
                from .extraction import extract_dependency

                for mask in solver.dependencies[:index]:
                    trial = extract_dependency(
                        prepared,
                        mask,
                        budget=budget,
                        product_policy=getattr(
                            collector_config, "extraction_policy", "reference"
                        ),
                    )
                    if trial.divisor is not None:
                        raise ValueError(
                            "checkpoint skipped a successful dependency"
                        )
                    engine.extractor.trials.append(trial)

                engine.extractor.next_dependency = index

            job.stats["matrix_reconstructions"] += 1

    if dynamic_snapshot:
        # Decoded input and rebuild scratch cease to be job-owned on return.
        # Keep the original public policy/identity and release only this
        # transient reservation; partials and prepared solver state stay live.
        config.metadata_reserve - public_config.metadata_reserve
        job.config = public_config
        if job.engine is not None:
            job.engine.config = replace(
                job.engine.config,
            )
            job.engine.collector.config = job.engine.config
        config = public_config

    if job.active and job.engine is None:
        raise ValueError("active SIQS checkpoint has no collector")
    if (extended or dynamic_snapshot) and job.base is not None:
        from .capacity import capacity_report

        job.stats["capacity"] = capacity_report(job.base, config)
    if extended and job.divisor is None:
        reasons = {
            "coefficient_limit": config.coefficient_trials
            > saved_config.coefficient_trials,
            "families_exhausted": config.family_count
            > saved_config.family_count,
            "stalled_yield": config.max_stalled > saved_config.max_stalled,
            "trivial_dependency_limit": config.max_trivial
            > saved_config.max_trivial,
        }
        if reasons.get(job.finished_reason, False):
            job.finished_reason = None
            if job.engine is not None and reasons.get(
                job.engine.storage_reason, False
            ):
                job.engine.storage_reason = None
                job.engine.storage_solve_done = False
        if reasons.get(job.pending_final_reason, False):
            # The extended search can collect again. Any retained downstream
            # prefix still runs before the engine admits another relation.
            job.pending_final_reason = None
    # Setup rebuilds the factor base and may reset the result field. A saved
    # completed split was independently checked above and must survive it.
    if payload["divisor"] is not None:
        job.divisor = payload["divisor"]
    return job
