"""Source-informed C1 diagnostics; no production relation admission changes."""

import argparse
import hashlib
import json
import platform
import random
import subprocess
import time
from collections import Counter
from dataclasses import asdict
from math import log1p, prod
from pathlib import Path
from unittest.mock import patch

from ....common import utils
from ....execution.budget import Budget, BudgetExhaustedError
from ....qs import SIQSJob
from ....qs.linear_algebra import filter_matrix, matrix_workspace
from ....qs.sieve_collector import SieveCollector
from ...infrastructure.performance.performance_audit import fingerprint
from ...support.paths import (
    source_path,
)
from ..b1.b1_calibration import common_config, load_corpus, run_one
from ..phase_three.phase_three_reference import _rss_bytes
from .c1_feasibility import (
    HERE,
    ROOT,
    ResidualAudit,
    kernel,
    load_control,
    machine_window,
    recover,
    split_residual,
    validate_control,
    verify_square,
    xor_selected,
)

CONTROL = HERE / "inputs/controls/c1_followup_observation.json"
MEMORY = 512 * 2**20
PREFIXES = (512, 2048, 8192, 32768, 131072)


def load_followup():
    control = json.loads(source_path(CONTROL).read_text())
    for name, expected in control["files"].items():
        if (
            hashlib.sha256((source_path(HERE / name)).read_bytes()).hexdigest()
            != expected
        ):
            raise ValueError("changed follow-up input: " + name)
    _, originals, configs = load_control()
    corpus = load_corpus(
        HERE / "inputs/corpora/phase_three_p34_large_v3_corpus.json"
    )
    fixtures = []
    for digits in (40, 50, 60):
        cases = [
            f
            for f in corpus["fixtures"]
            if f["kind"] == "balanced"
            and f["digits"] == digits
            and f["split"] == "training"
        ]
        selected = (
            cases[:2]
            if digits == 50
            else [
                next(f for f in originals if f["digits"] == digits),
                cases[0],
            ]
        )
        fixtures.extend(dict(f, c1_index=i) for i, f in enumerate(selected))
    configs[50] = common_config(30000, 65536, 6, rows=8192, partials=8192)
    return control, fixtures, configs


def graph_pairs(records):
    """Use vertex 1 for full and SLP rows; retain repeated endpoints."""
    pairs = []
    for record in records:
        endpoints = record["lp"]
        if len(endpoints) > 2:
            raise ValueError("not a two-endpoint residual")
        pairs.append(tuple([1] * (2 - len(endpoints)) + list(endpoints)))
    return pairs


def fundamental_cycles(pairs, budget, *, memory_bytes=MEMORY):
    """Offline complete multigraph basis, including disconnected components.

    Construct a spanning forest first. Each nonforest edge plus its unique
    forest path is one independent cycle: its nonforest edge appears in no
    other basis member. The count E-V+C proves the basis spans the kernel.
    """
    count = len(pairs)
    if count > 65536:
        raise ValueError("diagnostic edge cap exceeded")
    # At most 2E vertices: five maps, adjacency lists and parent tuples,
    # plus two adjacency entries per forest edge. 4 KiB/edge covers these
    # object/reference stores; the input records have a separate reservation.
    reserve = 32768 + 4096 * count
    if reserve > memory_bytes:
        raise MemoryError("offline graph reservation")
    parents, sizes, adjacency, extra = {}, {}, {}, []

    def find(vertex):
        while vertex != parents[vertex]:
            budget.consume(1)
            vertex = parents[vertex]
        return vertex

    for index, (left, right) in enumerate(pairs):
        budget.consume(1)
        for vertex in (left, right):
            if vertex not in parents:
                parents[vertex], sizes[vertex] = vertex, 1
                adjacency[vertex] = []
        a, b = find(left), find(right)
        if a == b:
            extra.append(index)
            continue
        if sizes[a] < sizes[b]:
            a, b = b, a
        parents[b] = a
        sizes[a] += sizes[b]
        adjacency[left].append((right, index))
        adjacency[right].append((left, index))

    components = [v for v in parents if parents[v] == v]
    assert len(extra) == count - len(parents) + len(components)
    # Each cycle retains one immutable E-bit integer. Allow twice its
    # payload plus object overhead, and eight simultaneous temporary masks.
    mask_bytes = (count + 7) // 8
    reserve += len(extra) * (256 + 2 * mask_bytes) + 8 * mask_bytes
    if reserve > memory_bytes:
        raise MemoryError("offline cycle provenance reservation")
    previous, depth = {}, {}
    for root in components:
        previous[root], depth[root] = (root, None), 0
        stack = [root]
        while stack:
            vertex = stack.pop()
            budget.consume(1)
            for child, index in adjacency[vertex]:
                if child in previous:
                    continue
                previous[child], depth[child] = (
                    (vertex, index),
                    depth[vertex] + 1,
                )
                stack.append(child)

    cycles, covered, lengths, disconnected = [], 0, Counter(), 0
    slp_root = find(1) if 1 in parents else None
    for index in extra:
        left, right = pairs[index]
        disconnected += find(left) != slp_root
        mask = 1 << index
        while left != right:
            budget.consume(1 + (count + 63) // 64)
            if depth[left] < depth[right]:
                left, right = right, left
            left, edge = previous[left]
            mask ^= 1 << edge
        cycles.append(mask)
        covered |= mask
        lengths[mask.bit_count()] += 1
    degrees = Counter(
        endpoint for pair in pairs for endpoint in pair if endpoint != 1
    )
    return tuple(cycles), dict(
        vertices=len(parents),
        components=len(components),
        largest_component=max((sizes[v] for v in components), default=0),
        forest_edges=count - len(extra),
        cycle_rank=len(extra),
        cycles_disconnected_from_slp=disconnected,
        cycle_lengths=dict(lengths),
        edges_outside_cycles=count - covered.bit_count(),
        degree_one_vertices=sum(degree == 1 for degree in degrees.values()),
        max_prime_degree=max(degrees.values(), default=0),
        reserve=reserve,
    )


def check_cycle(pairs, mask):
    """Independently count incidence, without relying on forest state."""
    odd = set()
    while mask:
        bit = mask & -mask
        for endpoint in pairs[bit.bit_length() - 1]:
            if endpoint in odd:
                odd.remove(endpoint)
            else:
                odd.add(endpoint)
        mask ^= bit
    return not odd


def validate_records(records, n, budget):
    """Check exact atoms and certification independently before analysis."""
    if len(records) > 65536:
        raise ValueError("diagnostic atom cap exceeded")
    seen, proven = set(), set()
    for record in records:
        for name in ("u", "square", "residual"):
            value = record[name]
            if not isinstance(value, int) or abs(value).bit_length() > 8192:
                raise ValueError("unbounded diagnostic integer")
        norm_bits = (record["u"] ** 2 - n).bit_length()
        previous = 1
        for prime, exponent in record["exponents"]:
            if not (
                isinstance(prime, int)
                and isinstance(exponent, int)
                and previous < prime <= 1000000
                and exponent >= 1
                and exponent * (prime.bit_length() - 1) <= norm_bits
            ):
                raise ValueError("unbounded or malformed diagnostic exponents")
            previous = prime
        if len(record["lp"]) > 2 or any(
            not isinstance(p, int) or not 2 <= p <= 10**12
            for p in record["lp"]
        ):
            raise ValueError("unbounded diagnostic endpoints")
        budget.consume(1 + len(record["exponents"]) * n.bit_length())
        if "polynomial" in record:
            key = record["polynomial"], record["position"]
            if key in seen:
                raise ValueError("duplicate exact atomic identity")
            seen.add(key)
        if not (
            record["sign"] in (-1, 1)
            and record["square"] >= 1
            and prod(record["lp"]) == record["residual"]
            and record["u"] ** 2 - n
            == record["sign"]
            * record["square"] ** 2
            * record["residual"]
            * prod(p**e for p, e in record["exponents"])
        ):
            raise ValueError("invalid exact diagnostic atom")
        for prime in record["lp"]:
            if prime not in proven:
                budget.consume(prime.bit_length() ** 2)
                if utils.classify_prime(prime) is not utils.Primality.PROVEN:
                    raise ValueError("uncertified diagnostic endpoint")
                proven.add(prime)


class AnalysisBudget:
    """Charge the parent ledger and impose a separate per-analysis deadline."""

    def __init__(self, parent):
        self.parent = parent
        self.local = Budget(work_limit=10**12, seconds=30, cpu_seconds=30)

    @property
    def used(self):
        return self.local.used

    def consume(self, amount=1):
        # Refusal may conservatively charge the parent, but never leaves
        # completed analysis work uncharged to the whole diagnostic run.
        if self.parent is not None:
            self.parent.consume(amount)
        self.local.consume(amount)


def compact_report(
    records, n, *, budget=None, memory_bytes=MEMORY, first_factor=False
):
    """Eliminate LPs sparsely and verify original-row square congruences."""
    budget = AnalysisBudget(budget)
    start_cpu, start_wall = time.process_time(), time.monotonic()
    try:
        validate_records(records, n, budget)
        pairs = graph_pairs(records)
        cycles, graph = fundamental_cycles(
            pairs, budget, memory_bytes=memory_bytes
        )
        base_labels = sorted({p for r in records for p, _ in r["exponents"]})
        columns = {p: i + 1 for i, p in enumerate(base_labels)}
        # The forest, all original parities, cycle masks, filter and solver
        # copies coexist. Retain the existing worst-case matrix reservation.
        reserve = graph["reserve"] + len(records) * (
            256 + 2 * ((len(columns) + 8) // 8)
        )
        reserve += 2 * matrix_workspace(len(cycles), len(columns) + 1)
        if reserve > memory_bytes:
            raise MemoryError("offline factor-base/provenance reservation")
        parity = []
        for record in records:
            row = int(record["sign"] < 0)
            for prime, exponent in record["exponents"]:
                if exponent % 2:
                    row ^= 1 << columns[prime]
            parity.append(row)
        rows = []
        for mask in cycles:
            budget.consume(mask.bit_count() * (1 + (len(records) + 63) // 64))
            assert check_cycle(pairs, mask)
            rows.append(xor_selected(parity, mask))
        rows = tuple(rows)
        filtered = filter_matrix(
            rows,
            weight_two=True,
            budget=budget,
            memory_bytes=memory_bytes - graph["reserve"],
        )
        # Verify every filter discovery and each lifted solver dependency.
        dependencies = list(filtered.zero_dependencies)
        for mask in kernel(filtered.rows, budget):
            dependencies.append(xor_selected(filtered.masks, mask))
        # Algebraic validation covers every generated mask. Extraction may
        # stop at a factor only for the separately frozen cost-attribution
        # audit; unextracted masks are never reported as square congruences.
        for mask in dependencies:
            budget.consume(mask.bit_count() * (1 + (len(columns) + 63) // 64))
            assert mask and xor_selected(rows, mask) == 0
        divisors, nontrivial, extracted = set(), 0, 0
        for mask in dependencies:
            lifted = xor_selected(cycles, mask)
            assert lifted and xor_selected(parity, lifted) == 0
            assert check_cycle(pairs, lifted)
            budget.consume(lifted.bit_count() * (1 + n.bit_length() ** 2))
            recovered = verify_square(records, lifted, n)
            extracted += 1
            nontrivial += bool(recovered)
            divisors.update(recovered)
            if first_factor and recovered:
                divisor = min(recovered)
                children = (divisor, n // divisor)
                budget.consume(
                    sum(child.bit_length() ** 2 for child in children)
                )
                if all(
                    utils.classify_prime(child)
                    is not utils.Primality.COMPOSITE
                    for child in children
                ):
                    break
        factors, remaining, labels = [], [n], []
        if divisors:
            divisor = min(divisors)
            remaining = []
            for child in (divisor, n // divisor):
                budget.consume(child.bit_length() ** 2)
                label = utils.classify_prime(child)
                if label is utils.Primality.COMPOSITE:
                    remaining.append(child)
                else:
                    factors.append(child)
                    labels.append(label.value)
        assert prod(factors) * prod(remaining) == n
        return dict(
            records=len(records),
            graph=graph,
            reserve=reserve,
            independent_lp_constraints=len(cycles),
            post_filter=filtered.stats,
            dependencies=extracted,
            algebraic_dependencies=len(dependencies),
            unextracted_dependencies=len(dependencies) - extracted,
            nontrivial_dependencies=nontrivial,
            trivial_dependencies=extracted - nontrivial,
            proper_divisors=sorted(divisors),
            factors=factors,
            remaining=remaining,
            certainty=labels,
            complete=not remaining,
            work=budget.used,
            cpu_seconds=time.process_time() - start_cpu,
            seconds=time.monotonic() - start_wall,
        )
    except (MemoryError, BudgetExhaustedError) as error:
        return dict(
            censored=str(error),
            records=len(records),
            factors=[],
            remaining=[n],
            certainty=[],
            complete=False,
            work=budget.used,
            cpu_seconds=time.process_time() - start_cpu,
            seconds=time.monotonic() - start_wall,
        )


class FollowupAudit(ResidualAudit):
    """Tight-product census with uniform sampling of rejected positions."""

    def __init__(self, config, seed, seconds):
        super().__init__(config, seed)
        if not (
            config.collector.division == "bucket"
            and config.collector.score_policy == "powers"
            and config.collector.score_backend == "list"
            and config.collector.small_prime_cutoff == 0
            and config.collector.threshold_extra == 0
        ):
            raise ValueError(
                "follow-up audit requires the frozen score configuration"
            )
        self.product_limit = 128 * config.base_bound**2
        self.inner_limit = 64 * config.base_bound**2
        self.budget = Budget(
            work_limit=10**13, seconds=seconds, cpu_seconds=seconds
        )
        self.sampling_rng = random.Random(548128 + seed)
        self.next_sample = self.sample_gap()
        self.sample_active = True
        self.costs = Counter()

    def sample_gap(self):
        # A geometric skip gives each position the same inclusion probability,
        # independent of block length. Floating arithmetic affects sampling
        # only; all factoring identities and threshold bounds remain integer.
        return int(log1p(-self.sampling_rng.random()) / log1p(-1 / 256))

    def classify(self, residual):
        if residual == 1:
            return (), "full"
        if residual > self.product_limit:
            return (), "product_bound"
        started = time.process_time()
        self.budget.consume(residual.bit_length() ** 2)
        prime = utils.classify_prime(residual) is utils.Primality.PROVEN
        classify_cost = time.process_time() - started
        self.classify_cpu += classify_cost
        self.costs["classify_128"] += classify_cost
        if residual <= self.inner_limit:
            self.costs["classify_64"] += classify_cost
        if prime:
            return (
                ((residual,), "slp")
                if residual <= (self.config.collector.residual_bound)
                else ((), "prime_above_slp")
            )
        if self.counts["split_attempts"] >= 131072:
            self.stop = "split_attempt_limit"
            return (), "split_unexamined"
        self.counts["split_attempts"] += 1
        started = time.process_time()
        endpoints, reason, evaluations = split_residual(
            residual, self.prime_limit, self.budget, seed=self.seed
        )
        cost = time.process_time() - started
        self.split_cpu += cost
        self.costs["split_128"] += cost
        if residual <= self.inner_limit:
            self.costs["split_64"] += cost
        self.counts["rho_evaluations"] += evaluations
        return endpoints, reason

    def _observe(self, collector, lo, hi, threshold):
        if self.stop:
            return
        self.budget.consume(hi - lo)
        self.counts["blocks"] += 1
        lower, _ = collector._bounds(lo, hi)
        wide_threshold = max(
            0, lower.bit_length() - 1 - (self.product_limit - 1).bit_length()
        )
        sampled = set()
        while self.sample_active and self.next_sample < hi - lo:
            sampled.add(self.next_sample)
            self.next_sample += 1 + self.sample_gap()
            if self.counts["uniform_positions"] + len(sampled) >= 32768:
                self.sample_active = False
                self.counts["sampling_prefix_positions"] = (
                    self.counts["positions"] + max(sampled) + 1
                )
        self.next_sample -= hi - lo
        for offset in range(hi - lo):
            self.budget.consume(1)
            self.counts["positions"] += 1
            score = collector._scores[offset]
            slp_block = score >= threshold
            self.counts[
                "slp_block_pass" if slp_block else "slp_block_reject"
            ] += 1
            wide = score >= wide_threshold
            if not wide and offset not in sampled:
                continue
            position = lo + offset
            value = collector.polynomial.value(position)
            slp_refine = collector._candidate_passes(value, offset)
            wide_refine = score >= max(
                0,
                abs(value).bit_length()
                - 1
                - (self.product_limit - 1).bit_length(),
            )
            if wide:
                self.counts[
                    "wide_refine_pass" if wide_refine else "wide_refine_reject"
                ] += 1
            if not (wide and wide_refine) and offset not in sampled:
                continue
            residual, exponents = recover(collector, position, offset)
            self.budget.consume((len(exponents) + 1) * abs(value).bit_length())
            if offset in sampled:
                self.budget.consume(
                    len(collector.factor_base.entries)
                    * max(1, abs(value).bit_length())
                )
                assert (residual, exponents) == recover(
                    collector, position, offset, full=True
                )
                self.counts["uniform_positions"] += 1
                sample_class = (
                    "uniform_pass"
                    if slp_block and slp_refine
                    else "uniform_reject"
                )
                self.sample(sample_class, residual)
                self.bits[f"{sample_class}:{residual.bit_length()}"] += 1
                if residual and residual <= self.product_limit:
                    assert wide and wide_refine
            if not (wide and wide_refine):
                continue
            self.bits[f"census:{residual.bit_length()}"] += 1
            if not residual:
                self.counts["zero"] += 1
                continue
            endpoints, reason = self.classify(residual)
            self.counts[reason] += 1
            if reason not in ("full", "slp", "dlp"):
                self.sample(reason, residual)
                if self.stop:
                    return
                continue
            if reason in ("full", "slp"):
                assert slp_block and slp_refine
            else:
                stage = (
                    "block"
                    if not slp_block
                    else (
                        "refinement" if not slp_refine else "residual_policy"
                    )
                )
                self.counts["dlp_lost_at_slp_" + stage] += 1
            reserve = 2048 + 128 * len(exponents)
            if (
                len(self.records) >= 65536
                or self.retained_bytes + reserve > 256 * 2**20
            ):
                self.stop = "diagnostic_record_cap"
                return
            polynomial = collector.polynomial
            self.records.append(
                dict(
                    kind=reason,
                    polynomial=polynomial.identity,
                    position=position,
                    u=polynomial.u_value(position),
                    sign=-1 if value < 0 else 1,
                    square=polynomial.square_coefficient,
                    exponents=exponents,
                    residual=residual,
                    lp=endpoints,
                    block=self.counts["blocks"],
                )
            )
            self.retained_bytes += reserve

    def snapshot(self, n):
        reports = {}
        for policy in ("slp", "64", "128"):
            records = [
                r
                for r in self.records
                if r["kind"] != "dlp"
                or (
                    policy != "slp"
                    and r["residual"]
                    <= int(policy) * self.config.base_bound**2
                )
            ]
            reports[policy] = compact_report(records, n, budget=self.budget)
        return dict(
            blocks=self.counts["blocks"],
            counts=dict(self.counts),
            costs=dict(self.costs),
            reports=reports,
            stopped=self.stop,
            retained_bytes=self.retained_bytes,
            diagnostic_work=self.budget.used,
        )


def save_capture(path, value):
    data = json.dumps(value, separators=(",", ":")) + "\n"
    if len(data.encode()) > 64 * 2**20:
        raise ValueError("follow-up capture cap exceeded")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        stream.write(data)


def followup_probe(fixture, seed, config, seconds):
    audit = FollowupAudit(config, seed, seconds)
    # Preserve three analysis slots plus bounded retained-prefix refinement.
    # Otherwise a terminal collection timeout would censor already retained
    # evidence without ever examining its final useful-yield opportunity.
    budget = Budget(
        work_limit=10**13,
        seconds=seconds - 120,
        cpu_seconds=seconds - 120,
        cancelled=lambda: audit.stop is not None,
    )
    job = SIQSJob(fixture["n"], seed=seed, config=config, budget=budget)
    original = SieveCollector._sieve

    def observed(collector, lo, hi, stats):
        threshold = original(collector, lo, hi, stats)
        audit.observe(collector, lo, hi, threshold)
        return threshold

    snapshots, prior = [], 0
    with patch.object(SieveCollector, "_sieve", observed):
        for blocks in PREFIXES:
            result = job.run(max_blocks=blocks - prior)
            prior = blocks
            snapshots.append(audit.snapshot(fixture["n"]))
            print(
                fixture["id"],
                audit.counts["blocks"],
                {
                    k: (
                        v.get("dependencies"),
                        v.get("proper_divisors"),
                        v.get("censored"),
                    )
                    for k, v in snapshots[-1]["reports"].items()
                },
                flush=True,
            )
            if result.reason != "paused" or audit.stop:
                break
    refine_prefixes(audit, snapshots, fixture["n"])
    assert (result.divisor or 1) * result.cofactor == fixture["n"]
    return dict(
        id=fixture["id"],
        digits=fixture["digits"],
        n=fixture["n"],
        seed=seed,
        config=asdict(config),
        snapshots=snapshots,
        rejection_samples=audit.samples,
        sample_populations=dict(audit.sample_seen),
        residual_bits=dict(audit.bits),
        sample_active=audit.sample_active,
        counts=dict(audit.counts),
        stopped=audit.stop,
        slp_result=dict(
            divisor=result.divisor,
            cofactor=result.cofactor,
            reason=result.reason,
            stats=result.stats,
        ),
        instrumented_seconds=budget.wall_used,
        instrumented_cpu=budget.cpu_used,
        production_work=budget.used,
        rss_bytes=_rss_bytes(),
    ), audit.records


class RefinementBudget:
    """A small suballowance, charged to the whole diagnostic run."""

    def __init__(self, parent):
        self.parent = parent
        self.local = Budget(work_limit=10**11, seconds=30, cpu_seconds=30)

    def consume(self, amount=1):
        self.parent.consume(amount)
        self.local.consume(amount)


def refine_prefixes(audit, snapshots, n):
    """Seek a directly verified DLP lead without collecting more positions."""
    bracket = None
    for policy in ("128", "64"):
        previous = 0
        for prefix in snapshots:
            result = prefix["reports"][policy]
            if "censored" in result:
                break
            if result["complete"]:
                bracket = previous, prefix["blocks"], prefix, policy
                break
            previous = prefix["blocks"]
        if bracket:
            break
    if bracket is None:
        return
    lower, upper, enclosing, policy = bracket
    budget = RefinementBudget(audit.budget)
    for _ in range(8):
        if upper - lower <= 1:
            break
        midpoint = (lower + upper) // 2
        prefix_records = [r for r in audit.records if r["block"] <= midpoint]
        reports = {}
        for candidate in ("slp", "64", "128"):
            selected = [
                r
                for r in prefix_records
                if r["kind"] != "dlp"
                or (
                    candidate != "slp"
                    and r["residual"]
                    <= int(candidate) * audit.config.base_bound**2
                )
            ]
            reports[candidate] = compact_report(selected, n, budget=budget)
        snapshots.append(
            dict(
                blocks=midpoint,
                reports=reports,
                costs=enclosing["costs"],
                stopped=None,
                refined=True,
                enclosing_blocks=enclosing["blocks"],
            )
        )
        if any("censored" in result for result in reports.values()):
            break
        if reports[policy]["complete"]:
            upper = midpoint
        else:
            lower = midpoint


def cell_decision(report, baseline, fixture):
    decisions = {}
    for policy in ("64", "128"):
        qualifying = []
        cumulative_analysis = 0.0
        for prefix in report["snapshots"]:
            cumulative_analysis += sum(
                r["cpu_seconds"] for r in prefix["reports"].values()
            )
            candidate, slp = (
                prefix["reports"][policy],
                prefix["reports"]["slp"],
            )
            for result in prefix["reports"].values():
                for divisor in result.get("proper_divisors", ()):
                    assert utils.valid_divisor(divisor, fixture["n"])
                    assert (
                        sorted((divisor, fixture["n"] // divisor))
                        == fixture["factors"]
                    )
            if (
                "censored" in candidate
                or "censored" in slp
                or prefix["stopped"]
            ):
                continue
            costs = prefix["costs"]
            charged = costs.get("split_" + policy, 0) + costs.get(
                "classify_" + policy, 0
            )
            charged += cumulative_analysis
            if (
                candidate["complete"]
                and candidate["proper_divisors"]
                and candidate["independent_lp_constraints"]
                - slp["independent_lp_constraints"]
                >= 16
                and charged
                <= max(baseline["cap_seconds"], baseline["cpu_seconds"]) / 2
            ):
                qualifying.append(
                    dict(
                        blocks=prefix["blocks"],
                        slp_has_factor=bool(slp["proper_divisors"]),
                        charged_cpu=charged,
                    )
                )
        decisions[policy] = qualifying
    return decisions


def campaign_decision(cells):
    """Require both independent inputs and an observed SLP-only shortfall."""
    indexed = {cell["name"]: cell["decision"] for cell in cells}
    passed = {}
    for digits in (40, 50, 60):
        policies = []
        for policy in ("64", "128"):
            witnesses = [
                indexed.get(f"{digits}-{i}", {}).get(policy, [])
                for i in (0, 1)
            ]
            if all(witnesses) and any(
                not prefix["slp_has_factor"]
                for case in witnesses
                for prefix in case
            ):
                policies.append(policy)
        if policies:
            passed[str(digits)] = policies
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--retained", type=Path)
    args = parser.parse_args()
    if (
        platform.python_implementation() != "PyPy"
        or platform.python_version_tuple()[:2] != ("3", "11")
    ):
        raise RuntimeError("C1 requires PyPy implementing Python 3.11")
    with machine_window():
        started, cpu = time.monotonic(), time.process_time()
        control, fixtures, configs = load_followup()
        summary = dict(
            control_sha256=hashlib.sha256(
                source_path(CONTROL).read_bytes()
            ).hexdigest(),
            sources=fingerprint(ROOT),
            driver_sha256=hashlib.sha256(
                source_path(Path(__file__)).read_bytes()
            ).hexdigest(),
            commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            runtime=platform.python_version(),
            cells=[],
        )
        save_capture(args.output / "manifest.json", summary)
        if args.retained:
            results = []
            for path in sorted(args.retained.glob("*-records.json")):
                original = json.loads(
                    source_path(
                        path.with_name(path.name.replace("-records", "-probe"))
                    ).read_text()
                )
                records = json.loads(source_path(path).read_text())
                cell_budget = Budget(
                    work_limit=10**12, seconds=30, cpu_seconds=30
                )
                results.append(
                    dict(
                        name=path.name,
                        source_sha256=hashlib.sha256(
                            source_path(path).read_bytes()
                        ).hexdigest(),
                        slp=compact_report(
                            [r for r in records if r["kind"] != "dlp"],
                            original["n"],
                            budget=cell_budget,
                        ),
                        dlp=compact_report(
                            records, original["n"], budget=cell_budget
                        ),
                    )
                )
            save_capture(args.output / "retained.json", results)
            return
        for fixture in fixtures:
            digits, index = fixture["digits"], fixture["c1_index"]
            seconds = {40: 180, 50: 480, 60: 900}[digits]
            control_seconds = {40: 30, 50: 120, 60: 300}[digits]
            if (
                max(time.monotonic() - started, time.process_time() - cpu)
                + seconds
                + control_seconds
                > 4320
            ):
                summary["stopped"] = "campaign_limit"
                break
            seed, name = (7, 29)[index], f"{digits}-{index}"
            with patch(
                "v2.benchmarks.qs.b1.b1_calibration.validate", validate_control
            ):
                baseline = run_one(
                    fixture, seed, configs[digits], control_seconds
                )
            save_capture(args.output / (name + "-control.json"), baseline)
            report, records = followup_probe(
                fixture, seed, configs[digits], seconds
            )
            report["decision"] = cell_decision(report, baseline, fixture)
            save_capture(args.output / (name + "-probe.json"), report)
            save_capture(args.output / (name + "-records.json"), records)
            summary["cells"].append(
                dict(name=name, decision=report["decision"])
            )
            if (
                sum(p.stat().st_size for p in args.output.glob("*.json"))
                > 512 * 2**20
            ):
                summary["stopped"] = "campaign_capture_limit"
                break
        summary["go_policies"] = campaign_decision(summary["cells"])
        summary.update(
            seconds=time.monotonic() - started,
            cpu_seconds=time.process_time() - cpu,
        )
        save_capture(args.output / "decision.json", summary)


if __name__ == "__main__":
    main()
