"""Experimental coarse SIQS collection with parent-owned finite allowances."""

import hashlib
import json
import multiprocessing
import os
import resource
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait
from dataclasses import asdict, dataclass, field, replace
from types import SimpleNamespace

from .. import arithmetic, utils
from ..arithmetic import pow
from ..budget import Budget, BudgetExhaustedError
from ..work_budget import PollingBudget
from .checkpoint import _restore_store, _store
from .factor_base import FactorBase, build_factor_base, checked_target
from .families import (
    PolynomialFamily,
    _checked_resources,
    _checksum,
    family_assignments,
)
from .pipeline import QSJob, QSResult
from .polynomial import Polynomial, qs_polynomial
from .relations import AtomicRelation, verify_atomic
from .sieve_collector import SieveCollector, SieveConfig

_PROCESS_STATE = None


@dataclass(frozen=True)
class ParallelConfig:
    """One fixed SIQS family schedule; no width recovery or auto dispatch.

    Memory includes parent, all worker workspaces, three bounded copies of
    each returned batch, and checkpoint encoding. Work leases are charged
    before submission; only reported unspent work is refunded.
    """

    backend: str = field(default="python-int", kw_only=True)
    base_bound: int = 200
    half_width: int = 256
    factor_count: int = 1
    family_count: int = 16
    pool_size: int = 32
    assignment_work: int = 10_000_000
    max_batch_atoms: int = 1024
    batch_width: int = 0
    poll_interval: int = 64
    parent_memory_bytes: int = 64 * 1024 * 1024
    worker_memory_bytes: int = 32 * 1024 * 1024
    memory_bytes: int = 512 * 1024 * 1024
    checkpoint_bytes: int = 4 * 1024 * 1024
    collector: SieveConfig = field(
        default_factory=lambda: SieveConfig(
            score_policy="powers",
            division="bucket",
            residual_bound=10000,
            max_atoms=8192,
            max_relations=4096,
            max_partials=1024,
        )
    )

    def __post_init__(self):
        if self.backend not in ("python-int", "gmpy2-mpz"):
            raise ValueError("unknown arithmetic backend")

        for name, low, high in (
            ("base_bound", 3, 100000),
            ("half_width", 1, 8192),
            ("factor_count", 1, 8),
            ("family_count", 1, 64),
            ("pool_size", self.factor_count, 128),
            ("assignment_work", 1, 2**63 - 1),
            ("max_batch_atoms", 1, 8192),
            ("batch_width", 0, 4096),
            ("poll_interval", 1, 64),
            ("parent_memory_bytes", 4 * 2**20, 2**30),
            ("worker_memory_bytes", 4 * 2**20, 2**30),
            ("memory_bytes", 4 * 2**20, 2**30),
            ("checkpoint_bytes", 4096, 16 * 2**20),
        ):
            value = utils.require_integer(getattr(self, name), name, low)
            if value > high:
                raise ValueError(f"{name} exceeds the parallel limit")

        if not isinstance(self.collector, SieveConfig):
            raise TypeError("collector must be a SieveConfig")
        if 8 * self.checkpoint_bytes >= self.parent_memory_bytes:
            raise MemoryError("checkpoint reservation exceeds parent memory")


def _rss():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _initialize_process(state):
    global _PROCESS_STATE
    with state["slot"].get_lock():
        index = state["slot"].value
        state["slot"].value += 1
    _PROCESS_STATE = (state, index)
    state["cpu"][index] = time.process_time()
    state["rss"][index] = _rss()


def _snapshot_process():
    """Publish post-transfer CPU; the barrier covers every live worker."""
    state, index = _PROCESS_STATE
    state["cpu"][index] = time.process_time()
    state["rss"][index] = _rss()
    state["barrier"].wait(timeout=30)


class _ParentBudget(PollingBudget):
    """Poll aggregate CPU during parent verification and matrix work."""

    def __init__(self, budget, pool, interval=64):
        super().__init__(budget, interval=interval, poll=pool.poll)


class _ExportCollector(SieveCollector):
    """Export verified atoms; matching belongs to the parent."""

    def _admit(self, atom, stats):
        if atom.relation_id in self._atoms:
            stats["duplicates"] += 1
            return None
        if len(self._atoms) >= self.config.max_atoms:
            return "batch_limit"
        reserve = 4096 + 256 * len(atom.exponents)
        reserve += 16 * (
            abs(atom.position).bit_length() + atom.polynomial.n.bit_length()
        )
        if self._workspace + reserve > self.config.memory_bytes:
            return "memory_limit"
        self.budget.consume(len(atom.exponents) + 1)
        self._workspace += reserve
        self._atoms[atom.relation_id] = atom
        stats["admitted_atoms"] += 1
        return None


class _FixedExportCollector(_ExportCollector):
    """Record direct splits while scanning the complete fixed workload."""

    found_divisor = None

    def _divide(self, position, offset, stats):
        atom, divisor = super()._divide(position, offset, stats)
        if divisor is not None:
            if self.found_divisor is None:
                self.found_divisor = divisor
            # Nonunit residuals cannot enter the central matching store.
            return None, None
        return atom, divisor


def _collect(task, state=None):
    """Return one complete family or an explicitly unpublished refusal."""
    base_record, primes, index, config, lease, deadline, fixed_work = task
    process = state is None
    slot = None
    if process:
        state, slot = _PROCESS_STATE
    started, cpu_started = time.perf_counter(), time.thread_time()
    last_poll = 0.0

    def cancelled():
        nonlocal last_poll
        if state["stop"].is_set():
            return True
        now = time.monotonic()
        if now - last_poll < 0.001:
            return False
        last_poll = now
        if process:
            state["cpu"][slot] = time.process_time()
        external = (
            sum(
                max(0, value - prior)
                for value, prior in zip(state["cpu"], state["cpu_base"])
            )
            if process
            else 0
        )
        if state["parent_cpu"].value + external >= state["cpu_limit"].value:
            state["stop"].set()
            return True
        return False

    budget = PollingBudget(
        Budget(
            work_limit=lease,
            seconds=None
            if deadline is None
            else max(0, deadline - time.monotonic()),
            cpu_seconds=None,
            cancelled=cancelled,
        ),
        interval=config.poll_interval,
    )
    collector, divisor, reason, scanned = None, None, "complete", 0
    peak = 0

    try:
        # Same-process jobs already own this checked immutable instance.
        # Serialized process records still require complete ingress checks.
        if isinstance(base_record, FactorBase):
            if process:
                raise ValueError("process workers require a serialized base")
            budget.consume(
                len(base_record.entries) * base_record.bound.bit_length() ** 2
            )
            base = base_record
        else:
            budget.consume(
                len(base_record[3]) * base_record[2].bit_length() ** 2
            )
            record = list(base_record)
            record[0] = arithmetic.get_backend(config.backend).integer(
                record[0]
            )
            base = FactorBase(*record)

        family = PolynomialFamily(
            base,
            primes,
            budget=budget,
            memory_bytes=config.worker_memory_bytes,
        )
        # The family and collector share this exact immutable base. Reserve
        # its storage once while their private metadata/buffers coexist.
        memory = (
            config.worker_memory_bytes
            - family.workspace_bytes
            + base.workspace_bytes
        )
        options = replace(
            config.collector,
            memory_bytes=max(0, memory),
            max_atoms=config.max_batch_atoms,
        )
        if config.batch_width:
            blocks = (
                2 * config.half_width + config.batch_width
            ) // config.batch_width
            gray_index, block = divmod(index % (family.count * blocks), blocks)
            steps = (family._make_step(gray_index, None),)
            lo = -config.half_width + block * config.batch_width
            hi = min(config.half_width + 1, lo + config.batch_width)
        else:
            steps = iter(family.next, None)
            lo, hi = -config.half_width, config.half_width + 1

        for step in steps:
            if collector is None:
                constructor = (
                    _FixedExportCollector if fixed_work else _ExportCollector
                )
                collector = constructor(
                    step.polynomial,
                    base,
                    config=options,
                    budget=budget,
                    precomputed_roots=step.roots,
                )
            else:
                collector.set_polynomial(step.polynomial, step.roots)

            result = collector.collect(lo, hi)
            scanned += result.stats["scanned"]
            peak = max(
                peak,
                result.workspace_bytes
                + family.workspace_bytes
                - base.workspace_bytes,
            )
            if result.reason != "complete":
                reason, divisor = result.reason, result.divisor
                break

        budget.consume(0)
    except BudgetExhaustedError:
        reason = budget.reason
    except MemoryError:
        reason = "memory_limit"
    finally:
        if process:
            state["cpu"][slot] = time.process_time()
            state["rss"][slot] = _rss()

    # Incomplete private prefixes are intentionally discarded and retried
    # under the same assignment ID. Their consumed work remains charged.
    atoms = (
        tuple(collector._atoms.values())
        if collector and reason == "complete"
        else ()
    )
    if fixed_work and reason == "complete":
        divisor = collector.found_divisor
    return dict(
        id=index,
        atoms=atoms,
        divisor=divisor,
        reason=reason,
        work=budget.used,
        scanned=scanned,
        workspace_bytes=peak,
        cpu_seconds=time.thread_time() - cpu_started,
        seconds=time.perf_counter() - started,
        pid=os.getpid(),
    )


class CollectionPool:
    """Reusable bounded executor; single active job, explicit close required.

    Spawned PyPy processes retain no relation stores between assignments.
    CPU publication occurs at most every 64 worker atomic actions (with a
    1 ms throttle). Parent accounting also polls every 64 atomic actions.
    Overshoot includes those chunks, native operations and the 10 ms wait;
    this is not an OS CPU/RSS sandbox. No nested backend threads are started.
    """

    def __init__(self, mode="serial", workers=1):
        if mode not in ("serial", "thread", "process"):
            raise ValueError("unknown execution mode")
        utils.require_integer(workers, "workers", 1)
        if workers > 4 or (mode == "serial" and workers != 1):
            raise ValueError("use one serial or at most four workers")
        if sys.implementation.name != "pypy" or sys.version_info[:2] != (
            3,
            11,
        ):
            raise RuntimeError("parallel collection requires PyPy Python 3.11")

        self.mode, self.workers = mode, workers
        context = None
        if mode == "process":
            context = multiprocessing.get_context("spawn")
            self.state = dict(
                stop=context.Event(),
                cpu=context.Array("d", workers),
                rss=context.Array("q", workers),
                cpu_base=context.Array("d", workers),
                parent_cpu=context.Value("d", 0),
                cpu_limit=context.Value("d", float("inf")),
                slot=context.Value("i", 0),
                barrier=context.Barrier(workers),
            )
        else:
            self.state = dict(
                stop=threading.Event(),
                cpu=[0.0] * workers,
                rss=[0] * workers,
                cpu_base=[0.0] * workers,
                parent_cpu=SimpleNamespace(value=0.0),
                cpu_limit=SimpleNamespace(value=float("inf")),
                slot=SimpleNamespace(value=0),
                barrier=None,
            )

        self.executor = None
        if mode == "process":
            self.executor = ProcessPoolExecutor(
                max_workers=workers,
                mp_context=context,
                initializer=_initialize_process,
                initargs=(self.state,),
            )
        elif mode == "thread":
            self.executor = ThreadPoolExecutor(max_workers=workers)

        self.lock, self.closed = threading.Lock(), False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self.state["stop"].set()
        if self.executor is not None:
            self.executor.shutdown(wait=True, cancel_futures=True)
        self.closed = True

    def begin(self, budget):
        if self.closed or not self.lock.acquire(blocking=False):
            raise ValueError("pool is closed or already running a job")
        self.state["stop"].clear()
        for i, value in enumerate(self.state["cpu"]):
            self.state["cpu_base"][i] = value
        self.state["cpu_limit"].value = (
            float("inf") if budget.cpu_seconds is None else budget.cpu_seconds
        )
        self.accounted_cpu = 0
        try:
            self.poll(budget)
        except BaseException:
            self.lock.release()
            raise

    def finish(self, budget):
        if self.mode == "process" and self.state["slot"].value:
            futures = [
                self.executor.submit(_snapshot_process)
                for _ in range(self.workers)
            ]
            for future in futures:
                future.result()

        self.poll(budget)

    def poll(self, budget):
        """Charge all published worker CPU, including spawn/import startup."""
        if self.mode == "process":
            current = sum(
                max(0, x - y)
                for x, y in zip(self.state["cpu"], self.state["cpu_base"])
            )
            budget.prior_cpu += max(0, current - self.accounted_cpu)
            self.accounted_cpu = current

        self.state["parent_cpu"].value = budget.cpu_used - self.accounted_cpu
        budget.consume(0)


class ParallelSIQSJob:
    """Ordered family merge and checked extraction, independent of dispatch.

    run drains/cancels all workers before returning. Completed pending batches
    and their admission cursor survive pauses; incomplete families replay from
    their beginning with the same ID and charged expenditure. Checkpoints do
    not store executors or raw positions. Solver caches rebuild on restoration.
    """

    def __init__(self, n, *, seed=7, config=None, budget=None):
        checked_target(n, 1)
        utils.require_integer(seed, "seed", 0)
        if seed >= 2**64:
            raise ValueError("seed exceeds 64 bits")
        self.n, self.seed = n, seed
        self.config = config if config is not None else ParallelConfig()
        if not isinstance(self.config, ParallelConfig):
            raise TypeError("config must be a ParallelConfig")
        self.n = arithmetic.get_backend(self.config.backend).integer(n)
        self.budget = (
            budget if budget is not None else Budget(work_limit=200_000_000)
        )
        self._original_budget = self.budget
        self._original_config = self._checked_config = self.config
        self.base = self.assignments = self.engine = None
        self.next_assignment, self.cursor = 0, 0
        self.pending = {}
        self._assignment_specs = {}
        self.refusals = {}
        self.divisor = None
        self.deferred_divisor = None
        self._running = False
        self.stats = dict(
            assignments=0,
            attempts=0,
            worker_work=0,
            cancelled_work=0,
            scanned=0,
            peak_owned_bytes=0,
            startup_and_collection_seconds=0.0,
            drain_seconds=0.0,
            worker_cpu_seconds=0.0,
        )

    def _setup(self):
        memory = (
            self.config.parent_memory_bytes - 8 * self.config.checkpoint_bytes
        )
        if self.base is None:
            result = build_factor_base(
                self.n,
                bound=self.config.base_bound,
                budget=self.budget,
                memory_bytes=memory,
            )
            self.base, self.divisor = result.factor_base, result.divisor
            if self.divisor:
                return

        if self.assignments is None:
            self.assignments = family_assignments(
                self.base,
                self.config.half_width,
                factor_count=self.config.factor_count,
                family_count=self.config.family_count,
                pool_size=self.config.pool_size,
                seed=self.seed,
                budget=self.budget,
                memory_bytes=memory,
            )

        if self.engine is None:
            self.engine = QSJob(
                qs_polynomial(self.base),
                self.base,
                0,
                0,
                budget=self.budget,
                weight_two=True,
                config=replace(self.config.collector, memory_bytes=memory),
            )

    @property
    def chunks_per_family(self):
        """Stable Gray/block schedule size; zero width retains coarse tasks."""
        if not self.config.batch_width:
            return 1
        blocks = (
            2 * self.config.half_width + self.config.batch_width
        ) // self.config.batch_width
        return (1 << (self.config.factor_count - 1)) * blocks

    @property
    def assignment_count(self):
        return len(self.assignments) * self.chunks_per_family

    def _assignment_primes(self, index):
        return self.assignments[index // self.chunks_per_family]

    def _memory(self, workers):
        # Atom shape is bounded by the base, input, family coefficients and
        # finite interval, before dispatch/pickling can allocate a result.
        # |A*F(x)| = |(A*x+B)**2-n|. A < bound**factor_count,
        # |B| <= A/2, and |x| <= half_width bound its bit length. Distinct
        # base-prime factors cannot outnumber the smallest-prime product
        # that fits that bound. Dense-base metadata is not attainable here.
        norm_bits = 1 + max(
            self.n.bit_length(),
            2
            * (
                self.config.factor_count * self.config.base_bound.bit_length()
                + self.config.half_width.bit_length()
                + 1
            ),
        )
        self.budget.consume(len(self.base.entries) + norm_bits)
        product, support = 1, 0
        for entry in self.base.entries:
            product *= entry.prime
            if product.bit_length() > norm_bits:
                break
            support += 1

        atom_bytes = 4096 + 256 * support
        atom_bytes += 128 * (
            self.n.bit_length()
            + 2
            * self.config.factor_count
            * self.config.base_bound.bit_length()
            + 32
        )
        positions = self.config.batch_width or (
            (2 * self.config.half_width + 1)
            * (1 << (self.config.factor_count - 1))
        )
        batch = min(self.config.max_batch_atoms, positions) * atom_bytes
        owned = self.config.parent_memory_bytes + workers * (
            self.config.worker_memory_bytes
            + 3 * batch
            + self.base.workspace_bytes
        )
        owned += 65536 + 1024 * len(self.assignments)
        if owned > self.config.memory_bytes:
            raise MemoryError(
                "aggregate worker/result reservations exceed memory_bytes"
            )
        self.stats["peak_owned_bytes"] = max(
            self.stats["peak_owned_bytes"], owned
        )

    def _merge(self, result):
        if result["id"] != self.next_assignment:
            raise ValueError("batch assignment order mismatch")
        if result["divisor"] is not None:
            if not utils.valid_divisor(result["divisor"], self.n):
                raise ValueError("invalid worker divisor")
            if not self.fixed_work:
                self.divisor = result["divisor"]
                return
            if self.deferred_divisor is None:
                self.deferred_divisor = result["divisor"]

        atoms = result["atoms"]
        if len(atoms) > self.config.max_batch_atoms:
            raise ValueError("oversized worker batch")
        stats = dict.fromkeys(
            (
                "duplicates",
                "dropped_partials",
                "evictions",
                "admitted_atoms",
                "matches",
            ),
            0,
        )

        while self.cursor < len(atoms):
            atom = atoms[self.cursor]
            self._check_atom_assignment(atom, self.next_assignment)
            verify_atomic(
                atom,
                self.base,
                residual_bound=self.config.collector.residual_bound,
                budget=self.engine.budget,
            )
            refusal = self.engine.collector._admit(atom, stats)
            if refusal:
                raise MemoryError(refusal)
            self.cursor += 1

        # A pending solve is finished before additional rows are admitted.
        if not self.fixed_work:
            self.divisor = self.engine._solve(False)
        del self.pending[self.next_assignment]
        self.next_assignment += 1
        self.cursor = 0
        self.stats["assignments"] += 1

    def _assignment_spec(self, index):
        """Cache at most eight immutable assignment/window specifications."""
        primes = self._assignment_primes(index)
        key = index, primes, self.config.half_width, self.config.batch_width
        if key in self._assignment_specs:
            return self._assignment_specs[key]
        a = 1
        roots = []

        for prime in primes:
            a *= prime
            roots.append(
                next(
                    e.square_roots[0]
                    for e in self.base.entries
                    if e.prime == prime
                )
            )

        b = None
        lo, hi = -self.config.half_width, self.config.half_width + 1
        if self.config.batch_width:
            blocks = (
                2 * self.config.half_width + self.config.batch_width
            ) // self.config.batch_width
            gray_index, block = divmod(index % self.chunks_per_family, blocks)
            gray = gray_index ^ (gray_index >> 1)
            raw = 0
            for offset, (prime, root) in enumerate(zip(primes, roots)):
                quotient = a // prime
                term = quotient * pow(quotient, -1, prime) * root % a
                raw += -term if offset and gray & (1 << (offset - 1)) else term
            b = (raw + a // 2) % a - a // 2
            lo += block * self.config.batch_width
            hi = min(hi, lo + self.config.batch_width)

        if len(self._assignment_specs) == 8:
            del self._assignment_specs[next(iter(self._assignment_specs))]
        spec = a, roots[0], b, lo, hi
        self._assignment_specs[key] = spec
        return spec

    def _check_atom_assignment(self, atom, index):
        a, first_root, b, lo, hi = self._assignment_spec(index)
        primes = self._assignment_primes(index)
        polynomial = atom.polynomial
        if (
            polynomial.a != a
            or abs(polynomial.b) > a // 2
            or polynomial.b % primes[0] != first_root
            or not lo <= atom.position < hi
            or (b is not None and polynomial.b != b)
        ):
            raise ValueError("relation outside its assigned family/Gray/block")

    def _atom_chunk(self, atom, products):
        """Locate provenance in bounded metadata, without scanning history."""
        family_index = products.get(atom.polynomial.a)
        if family_index is None:
            raise ValueError("store includes an unknown family")
        primes = self.assignments[family_index]
        gray = 0

        for offset, prime in enumerate(primes[1:]):
            root = next(
                e.square_roots[0]
                for e in self.base.entries
                if e.prime == prime
            )
            residue = atom.polynomial.b % prime
            if residue == (-root) % prime:
                gray |= 1 << offset
            elif residue != root:
                raise ValueError("store includes an unknown Gray polynomial")

        index, shifted = gray, gray >> 1
        while shifted:
            index ^= shifted
            shifted >>= 1
        block = (
            atom.position + self.config.half_width
        ) // self.config.batch_width
        blocks = (
            2 * self.config.half_width + self.config.batch_width
        ) // self.config.batch_width
        assignment = (
            family_index * (1 << (self.config.factor_count - 1)) + index
        ) * blocks + block
        if not 0 <= block < blocks:
            raise ValueError("store includes an unknown block")
        self._check_atom_assignment(atom, assignment)
        return assignment

    def _check_config(self):
        """Keep assignment identity and existing storage reservations fixed."""
        if self.config is self._checked_config:
            return
        if not isinstance(self.config, ParallelConfig):
            raise TypeError("config must be a ParallelConfig")
        original = self._original_config
        normalized = replace(
            self.config,
            assignment_work=original.assignment_work,
            poll_interval=original.poll_interval,
            checkpoint_bytes=original.checkpoint_bytes,
        )
        if normalized != original or (
            self.config.checkpoint_bytes > original.checkpoint_bytes
        ):
            raise ValueError(
                "changed assignment/storage settings require a new job"
            )

        self._checked_config = self.config

    def run(self, *, pool=None, max_assignments=None, fixed_work=False):
        """Run ordered collection; extend the original Budget to resume.

        fixed_work collects the entire finite schedule before extraction for
        throughput experiments. Default extraction checks each family batch.
        The mode must stay the same across resumptions.
        """
        self._check_config()
        if self.budget is not self._original_budget:
            raise ValueError("extend the original Budget in place")
        if self._running:
            raise ValueError("job is already running")
        if max_assignments is not None:
            utils.require_integer(max_assignments, "max_assignments", 1)
        if type(fixed_work) is not bool:
            raise TypeError("fixed_work must be Boolean")
        if hasattr(self, "fixed_work") and self.fixed_work != fixed_work:
            raise ValueError("resume must retain fixed_work mode")
        self.fixed_work = fixed_work
        if self.divisor is not None:
            return self._result("factor_found")
        own_pool = pool is None
        pool = pool if pool is not None else CollectionPool()
        inflight, reason, merged = {}, "families_exhausted", 0
        begun = False
        self._running = True

        try:
            pool.begin(self.budget)
            begun = True
            self._setup()
            slots = max(
                pool.workers,
                max(self.pending, default=self.next_assignment)
                - self.next_assignment
                + 1,
            )
            if not self.divisor:
                self._memory(slots)
                parent_budget = _ParentBudget(
                    self.budget, pool, self.config.poll_interval
                )
                self.engine.budget = self.engine.collector.budget = (
                    parent_budget
                )

            while (
                not self.divisor
                and self.next_assignment < self.assignment_count
            ):
                pool.poll(self.budget)
                try:
                    if (
                        self.engine.solver is not None
                        or self.engine.extractor is not None
                    ):
                        self.divisor = self.engine._solve(False)
                        if self.divisor:
                            break

                    if self.next_assignment in self.pending:
                        self._merge(self.pending[self.next_assignment])
                        merged += 1
                        if (
                            max_assignments is not None
                            and merged >= max_assignments
                        ):
                            reason = "paused"
                            break

                        continue
                except BudgetExhaustedError:
                    if self.budget.reason != "work_limit" or not inflight:
                        raise

                    # Running leases can refund unused work. Keep admission
                    # and solver cursors, then receive before trying again.
                if not inflight:
                    stop = min(
                        self.assignment_count,
                        self.next_assignment + pool.workers,
                    )
                    if max_assignments is not None:
                        stop = min(
                            stop,
                            self.next_assignment + max_assignments - merged,
                        )

                    started = time.perf_counter()

                    for index in range(self.next_assignment, stop):
                        if index in self.pending:
                            continue
                        available = self.budget.work_limit - self.budget.used
                        if available == 0 and inflight:
                            break
                        lease = min(self.config.assignment_work, available)
                        if lease == 0:
                            self.budget.consume(1)
                        # Reserve before dispatch so simultaneous workers
                        # cannot each spend the same remaining allowance.
                        self.budget.consume(lease)
                        self.stats["attempts"] += 1
                        seconds = self.budget.seconds
                        deadline = (
                            None
                            if seconds is None
                            else time.monotonic()
                            + max(0, seconds - self.budget.wall_used)
                        )
                        task = (
                            (
                                self.base.n,
                                self.base.multiplier,
                                self.base.bound,
                                self.base.entries,
                            )
                            if pool.mode == "process"
                            else self.base,
                            self._assignment_primes(index),
                            index,
                            self.config,
                            lease,
                            deadline,
                            fixed_work,
                        )
                        if pool.executor is None:
                            result = _collect(task, pool.state)
                            self._receive(result, lease)
                        else:
                            try:
                                future = (
                                    pool.executor.submit(_collect, task)
                                    if pool.mode == "process"
                                    else pool.executor.submit(
                                        _collect, task, pool.state
                                    )
                                )
                            except BaseException:
                                self.budget.used -= lease
                                raise

                            inflight[future] = lease

                    self.stats["startup_and_collection_seconds"] += (
                        time.perf_counter() - started
                    )
                    if (
                        pool.executor is None
                        and self.next_assignment not in self.pending
                    ):
                        reason = result["reason"]
                        break

                if inflight:
                    started = time.perf_counter()
                    done, _ = wait(inflight, timeout=0.01)
                    self.stats["startup_and_collection_seconds"] += (
                        time.perf_counter() - started
                    )
                    for future in done:
                        lease = inflight.pop(future)
                        self._receive(future.result(), lease)
                    if (
                        not inflight
                        and self.next_assignment not in self.pending
                    ):
                        reason = self.refusals.get(
                            self.next_assignment, "cancelled"
                        )
                        break

            if (
                not self.divisor
                and self.next_assignment == self.assignment_count
            ):
                self.divisor = self.engine._solve(True)
                if self.divisor is None:
                    self.divisor = self.deferred_divisor
                reason = "families_exhausted"
        except BudgetExhaustedError:
            reason = self.budget.reason
        except MemoryError:
            reason = "memory_limit"
        finally:
            try:
                if begun:
                    started = time.perf_counter()
                    pool.state["stop"].set()
                    errors = []

                    for future, lease in list(inflight.items()):
                        if future.cancel():
                            self.budget.used -= lease
                        else:
                            try:
                                self._receive(future.result(), lease)
                            except Exception as error:
                                self.stats["cancelled_work"] += lease
                                errors.append(error)

                    self.stats["drain_seconds"] += (
                        time.perf_counter() - started
                    )
                    try:
                        pool.finish(self.budget)
                    except BudgetExhaustedError:
                        if not self.divisor:
                            reason = self.budget.reason
                    except Exception as error:
                        errors.append(error)

                    if errors:
                        raise errors[0]
            finally:
                if begun:
                    pool.lock.release()
                    if self.engine is not None:
                        self.engine.budget = self.engine.collector.budget = (
                            self.budget
                        )

                        # A paused solver/extractor must not pin the old pool
                        # through a parent polling callback after draining it.
                        for pending in (
                            self.engine.solver,
                            self.engine.extractor,
                        ):
                            if pending is not None:
                                pending.budget = self.budget

                try:
                    if own_pool:
                        pool.close()
                finally:
                    self._running = False

        return self._result("factor_found" if self.divisor else reason)

    def _receive(self, result, lease):
        work = utils.require_integer(result["work"], "worker work", 0)
        if work > lease:
            raise ValueError("worker exceeded its work lease")
        # Return only unused reservations; interrupted work remains charged.
        self.budget.used -= lease - work
        self.stats["worker_work"] += work
        self.stats["worker_cpu_seconds"] += result["cpu_seconds"]
        self.stats["scanned"] += result["scanned"]
        if result["reason"] in ("complete", "factor_found"):
            index = result["id"]
            if index in self.pending or index < self.next_assignment:
                raise ValueError("duplicate worker assignment")
            self.pending[index] = result
        else:
            self.refusals[result["id"]] = result["reason"]
            self.stats["cancelled_work"] += work

    def _result(self, reason):
        if self.divisor is not None and not utils.valid_divisor(
            self.divisor, self.n
        ):
            raise ArithmeticError("invalid parallel split")
        stats = dict(self.stats)
        stats.update(
            work_used=self.budget.used,
            cpu_used=self.budget.cpu_used,
            pending_assignments=sorted(self.pending),
            schedule_complete=self.assignments is not None
            and self.next_assignment == self.assignment_count,
            chunks_per_family=self.chunks_per_family,
        )
        if self.engine is not None:
            stats["engine"] = self.engine.stats
        return QSResult(
            self.divisor,
            self.n // self.divisor if self.divisor else self.n,
            reason,
            self.next_assignment,
            stats,
        )

    def checkpoint(self):
        """Encode quiescent provenance without executable pickle payloads."""
        self._check_config()
        if self._running:
            raise ValueError("checkpoint requires drained workers")

        def encode(atom):
            return [
                atom.polynomial.a,
                atom.polynomial.b,
                atom.position,
                atom.sign,
                atom.exponents,
                atom.residual,
            ]

        payload = dict(
            version=4,
            backend=arithmetic.get_backend(self.config.backend).identity,
            n=self.n,
            seed=self.seed,
            config=asdict(self.config),
            next_assignment=self.next_assignment,
            cursor=self.cursor,
            divisor=self.divisor,
            deferred_divisor=self.deferred_divisor,
            fixed_work=getattr(self, "fixed_work", False),
            stats=self.stats,
            store=_store(self.engine.collector) if self.engine else None,
            pending=[
                dict(
                    id=i,
                    atoms=[encode(a) for a in r["atoms"]],
                    divisor=r["divisor"],
                    reason=r["reason"],
                )
                for i, r in sorted(self.pending.items())
            ],
        )
        parts, size = [], 0

        for part in json.JSONEncoder(
            default=arithmetic.json_integer,
            sort_keys=True,
            separators=(",", ":"),
        ).iterencode(payload):
            size += len(part.encode())
            if size + 1024 > self.config.checkpoint_bytes:
                raise MemoryError(
                    "parallel checkpoint exceeds checkpoint_bytes"
                )
            parts.append(part)

        blob = "".join(parts)
        resources = dict(
            work_used=self.budget.used,
            wall_used=self.budget.wall_used,
            cpu_used=self.budget.cpu_used,
        )
        return dict(
            blob=blob,
            sha256=hashlib.sha256(blob.encode()).hexdigest(),
            resources=resources,
            resources_sha256=_checksum(resources),
        )

    @classmethod
    def from_checkpoint(cls, checkpoint, *, budget):
        """Reverify retained atoms under the cumulative allowance."""
        blob = checkpoint["blob"]
        if not isinstance(blob, str) or len(blob.encode()) > 16 * 2**20:
            raise ValueError("invalid parallel checkpoint size")
        if hashlib.sha256(blob.encode()).hexdigest() != checkpoint["sha256"]:
            raise ValueError("parallel checkpoint digest mismatch")
        resources = _checked_resources(checkpoint["resources"])
        if _checksum(resources) != checkpoint["resources_sha256"]:
            raise ValueError("parallel checkpoint resource digest mismatch")
        payload = json.loads(blob)
        if type(payload["version"]) is not int or payload["version"] not in (
            1,
            2,
            3,
            4,
        ):
            raise ValueError("unknown parallel checkpoint version")

        if payload["version"] >= 3 and payload["store"] is not None:
            if not isinstance(payload["store"], dict) or (
                "row_order" not in payload["store"]
            ):
                raise ValueError("mixed-order checkpoint lacks row_order")

        options = payload["config"]
        options["collector"] = SieveConfig(**options["collector"])
        config = ParallelConfig(**options)
        if (
            payload.get(
                "backend", "python-int" if payload["version"] < 4 else None
            )
            != arithmetic.get_backend(config.backend).identity
        ):
            raise ValueError("incompatible checkpoint backend")
        if payload["version"] == 1 and config.batch_width:
            raise ValueError("legacy checkpoint cannot contain chunk progress")
        if len(blob.encode()) + 1024 > config.checkpoint_bytes:
            raise MemoryError("parallel checkpoint exceeds configured cap")
        budget.used = max(budget.used, resources["work_used"])
        budget.prior_wall = max(budget.prior_wall, resources["wall_used"])
        budget.prior_cpu = max(budget.prior_cpu, resources["cpu_used"])
        if budget.used > budget.work_limit:
            raise ValueError("work limit is below restored expenditure")
        job = cls(
            payload["n"], seed=payload["seed"], config=config, budget=budget
        )
        job._setup()
        next_index = utils.require_integer(
            payload["next_assignment"], "next assignment", 0
        )
        if job.assignments is None:
            if next_index or payload["store"] or payload["pending"]:
                raise ValueError("invalid setup-factor checkpoint")
        else:
            if next_index > job.assignment_count:
                raise ValueError("assignment progress exceeds schedule")
            if payload["store"] is not None:
                _restore_store(payload["store"], job.engine.collector, budget)
            pending = payload["pending"]
            if not isinstance(pending, list) or len(pending) > 4:
                raise ValueError("pending batch cap exceeded")
            job._memory(max(1, len(pending)))
            for record in pending:
                index = utils.require_integer(
                    record["id"], "assignment ID", next_index
                )
                if (
                    index >= min(job.assignment_count, next_index + 4)
                    or index in job.pending
                ):
                    raise ValueError("invalid pending assignment")

                if (
                    record["reason"] not in ("complete", "factor_found")
                    or len(record["atoms"]) > config.max_batch_atoms
                ):
                    raise ValueError("invalid pending batch")

                if record["divisor"] is not None and not utils.valid_divisor(
                    record["divisor"], job.n
                ):
                    raise ValueError("invalid pending checkpoint divisor")
                atoms = []

                for a, b, x, sign, exponents, residual in record["atoms"]:
                    atom = AtomicRelation(
                        Polynomial(job.n, 1, a, b),
                        x,
                        sign,
                        tuple(tuple(pair) for pair in exponents),
                        residual,
                    )
                    verify_atomic(
                        atom,
                        job.base,
                        residual_bound=config.collector.residual_bound,
                        budget=budget,
                    )
                    job._check_atom_assignment(atom, index)
                    atoms.append(atom)

                job.pending[index] = dict(record, atoms=tuple(atoms))

        job.next_assignment = next_index
        job.cursor = utils.require_integer(
            payload["cursor"], "admission cursor", 0
        )
        if job.cursor > len(job.pending.get(next_index, {}).get("atoms", ())):
            raise ValueError("cursor exceeds pending batch")
        if job.engine is not None:
            products = {}
            if config.batch_width:
                budget.consume(len(job.assignments) * config.factor_count)
                for index, primes in enumerate(job.assignments):
                    a = 1
                    for prime in primes:
                        a *= prime
                    products[a] = index

            prefix = {
                a.relation_id: a
                for a in job.pending.get(next_index, {}).get("atoms", ())[
                    : job.cursor
                ]
            }

            for atom in job.engine.collector._atoms.values():
                if config.batch_width:
                    index = job._atom_chunk(atom, products)
                    if index > next_index or (
                        index == next_index
                        and prefix.get(atom.relation_id) != atom
                    ):
                        raise ValueError("store includes an uncommitted chunk")

                    continue

                for index in range(next_index + bool(job.cursor)):
                    try:
                        job._check_atom_assignment(atom, index)
                        if (
                            index == next_index
                            and prefix.get(atom.relation_id) != atom
                        ):
                            raise ValueError(
                                "store includes an uncommitted atom"
                            )

                        break
                    except ValueError:
                        continue
                else:
                    raise ValueError("store includes an uncommitted family")

        job.divisor = payload["divisor"]
        job.deferred_divisor = payload["deferred_divisor"]
        if job.deferred_divisor is not None and not utils.valid_divisor(
            job.deferred_divisor, job.n
        ):
            raise ValueError("invalid deferred checkpoint divisor")
        if job.divisor is not None and not utils.valid_divisor(
            job.divisor, job.n
        ):
            raise ValueError("invalid checkpoint divisor")
        if type(payload["fixed_work"]) is not bool:
            raise ValueError("invalid fixed-work flag")
        job.fixed_work = payload["fixed_work"]
        if set(payload["stats"]) != set(job.stats):
            raise ValueError("invalid parallel statistics shape")
        for name, value in payload["stats"].items():
            if isinstance(job.stats[name], int):
                utils.require_integer(value, name, 0)
            elif (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not 0 <= value < float("inf")
            ):
                raise ValueError("invalid parallel timing")

        if payload["stats"]["assignments"] != next_index:
            raise ValueError("assignment progress disagrees with statistics")
        job.stats = payload["stats"]
        return job
