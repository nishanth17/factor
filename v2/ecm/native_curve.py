"""Whole ECM in ordinary GMP with canonical, durable scalar checkpoints."""

import random
import threading
from functools import wraps

from ..common import utils
from ..execution import stage_jobs
from .chain_records import verify_frontier
from .native_packet import NativePacket, _bytes
from .native_resources import close_resource_checks, poll_native_budget

PHASES = (
    "setup",
    "stage_one",
    "replay",
    "stage_two_setup",
    "baby_steps",
    "giant_setup",
    "stage_two",
    "term_replay",
)


def exclusive_curve(method):
    """Reject close or a second caller while C borrows this curve."""

    @wraps(method)
    def call(self, *args, **kwargs):
        if not self.call_lock.acquire(blocking=False):
            raise RuntimeError("native curve is active")
        try:
            return method(self, *args, **kwargs)
        finally:
            self.call_lock.release()

    return call


class BlockStore:
    """One checked immutable native metadata store per retained scalar/team."""

    def __init__(self):
        self.blocks = {}
        self.verified = set()
        self.used = 0
        self.ffi = self.lib = None

    def bind(self, packet):
        if self.ffi is None:
            self.ffi, self.lib = packet.ffi, packet.lib
        elif self.lib is not packet.lib:
            raise ValueError("native metadata belongs to another library")

    def close(self):
        close_resource_checks(self)
        for block, _ in self.blocks.values():
            if self.lib.b4_block_destroy(block):
                raise RuntimeError("native metadata still borrowed by a curve")
        self.blocks.clear()
        self.verified.clear()
        self.used = 0

    def register(self, key, values, lo, hi, bound, plan, budget, *, b1, b2):
        cached = self.blocks.get(key)
        if cached is not None:
            budget.consume(1)
            return cached
        values = tuple(values)
        if len(self.blocks) >= 4096 or len(values) > 65536:
            raise ValueError(
                "native prime block geometry exceeds supported limits"
            )
        reserve = 512 + 128 * len(values)
        if plan is not None:
            reserve += sum(
                2
                * (
                    len(plan.entries[p][1].record.code)
                    + len(plan.entries[p][2].record.code)
                )
                for p in values
            )
        if self.used + reserve > 64 * 2**20:
            raise ValueError(
                "native immutable metadata exceeds supported geometry"
            )
        budget.consume(reserve)
        ffi, lib = self.ffi, self.lib
        units = ffi.new("b4_curve_unit[]", max(1, len(values)))
        owners = [units]
        try:
            for target, prime in zip(units, values):
                target.prime = prime
                target.power = (
                    utils.prime_power(prime, bound) if bound else prime
                )
                if plan is None:
                    continue
                power, fast, unit, reserve_work = plan.entries[prime]
                if power != target.power:
                    raise ValueError(
                        "native program disagrees with prime power"
                    )
                target.reserve = reserve_work
                for destination, action in (
                    (target.fast, fast),
                    (target.unit, unit),
                ):
                    certificate = (action.record, action.masks)
                    if certificate not in self.verified:
                        verify_frontier(*certificate)
                        self.verified.add(certificate)
                    code = ffi.new("unsigned char[]", action.record.code)
                    masks = ffi.new("unsigned char[]", action.masks)
                    owners.extend((code, masks))
                    destination.code, destination.masks = code, masks
                    destination.code_size = len(action.record.code)
                    destination.masks_size = len(action.masks)
                    destination.slots, destination.output = (
                        action.record.slots,
                        action.record.output,
                    )
            block = lib.b4_block_create(
                lo, hi, bound or 0, int(plan is not None), units, len(values)
            )
            if block == ffi.NULL:
                raise ValueError(
                    "native immutable prime/program block refused"
                )
            try:
                if not bound and lib.b4_block_prepare(block, b1, b2):
                    raise ValueError(
                        "native prepared stage-two indices refused"
                    )
                self.blocks[key] = block, values
                self.used += reserve
            except BaseException:
                lib.b4_block_destroy(block)
                raise
        finally:
            for owner in reversed(owners):
                ffi.release(owner)
        # Registration uses the same exact root transport as native slices;
        # child trees and custom providers retain the original observation.
        poll_native_budget(budget, self)
        return self.blocks[key]


class NativeCurvePacket(NativePacket):
    """Private scratch with exports only at complete canonical actions."""

    full_ecm = True
    atomic_actions = False
    phase_clocks = False
    actions = None

    def __init__(self, *, mode="c-gmp", _blocks=None, **kwargs):
        super().__init__(mode=mode, **kwargs)
        self.blocks = _blocks or BlockStore()
        self.owns_blocks = _blocks is None
        self.job = self.context_source = self.config = None
        self.cursor_blocks = {}
        self.charged = 0
        self.calls = 0
        self.info = self.atom = self.pair = None
        self.plan = None
        self.statistics = dict(setup_ops=0, stage_one_ops=0, stage_two_ops=0)
        self.phase_times = [0, 0, 0]
        self.published = False
        self.call_lock = threading.Lock()
        self.bound_identity = None
        self.canonical_ready = False
        self.lifetime_statistics = dict(self.statistics)
        self.lifetime_phase_times = [0, 0, 0]

    def _buffers(self):
        if self.info is None:
            if not self.lib.b4_whole_ecm():
                raise ValueError("native library lacks whole-ECM capability")
            self.info = self.ffi.new("b4_curve_info *")
            self.atom = self.ffi.new("b4_atom *")
            self.pair = self.ffi.new("unsigned long long[2]")
            self.blocks.bind(self)

    def _read(self):
        if self.lib.b4_curve_info_get(self.context, self.info):
            raise ValueError("native curve has no checked state")
        return self.info[0]

    def _field(self, field, index=0):
        if self.lib.b4_curve_get(self.context, field, index, self.atom):
            raise ValueError("invalid native canonical field")
        size = int(self.atom.size)
        if not 1 <= size <= 128:
            raise ValueError("native curve returned an invalid integer length")
        raw = bytes(self.ffi.buffer(self.atom.data, size))
        if size > 1 and raw[0] == 0:
            raise ValueError("native curve returned a noncanonical integer")
        value = int.from_bytes(raw, "big")
        if value >= self.job["n"]:
            raise ValueError("native curve field exceeds its modulus")
        return value

    def _plan(self, job, context, budget):
        return stage_jobs.stage_one_plan(job, context, budget)

    def _block(self, lo, hi, bound, budget, *, values=None, purpose="cursor"):
        plan = self.plan if bound else None
        identity = None if plan is None else plan.identity
        key = (
            purpose,
            lo,
            hi,
            bound,
            identity,
            self.job["b1"],
            self.job["b2"],
        )
        if key not in self.blocks.blocks and values is None:
            program = getattr(self.context_source, "program_segment", None)
            if program is not None:
                values = program(lo, hi, budget, bound=bound)
            else:
                budget.consume(
                    self.context_source.segment_size
                    + len(self.context_source.base_primes)
                )
                values = list(self.context_source.primes(lo, hi))
        block = self.blocks.register(
            key,
            () if values is None else values,
            lo,
            hi,
            bound,
            plan,
            budget,
            b1=self.job["b1"],
            b2=self.job["b2"],
        )
        if purpose == "cursor":
            self.cursor_blocks[(lo, hi, bound)] = block
        return block

    def bind_job(self, job, budget, context, config):
        if self.closed:
            raise RuntimeError("native curve packet is closed")
        identity = (
            job["n"],
            job["seed"],
            job["b1"],
            job["b2"],
            config.chunk_size,
            config.gcd_batch,
            config.segment_size,
            config.ecm_chain_mode,
            config.ecm_chain_family,
        )
        if self.job is job:
            if self.bound_identity != identity:
                raise ValueError("live native curve identity changed")
            if self.canonical_ready:
                return
        if self.job is not None:
            self.discard(sync=not self.published and self.canonical_ready)
        if not self.available(job["n"], budget):
            raise ValueError("whole ECM requires an available native backend")
        self._buffers()
        if (
            config.ecm_pair_distance is not None
            or config.ecm_pair_wheel is not None
        ):
            raise ValueError(
                "whole native ECM supports classical continuation only"
            )
        if config.ecm_chain_family != "auto":
            raise ValueError(
                "whole native ECM requires the accepted PRAC/ladder policy"
            )
        self.job, self.context_source, self.config = job, context, config
        self.bound_identity = identity
        self.plan = None
        self.canonical_ready = False
        self.cursor_blocks.clear()
        self.charged = 0
        self.published = False
        sigma = random.Random(job["seed"]).randrange(6, 2**63)
        owners, values = self._values([job["n"], sigma])
        try:
            status = self.lib.b4_curve_begin(
                self.context,
                values[0][0],
                values[1][0],
                job["b1"],
                job["b2"],
                config.chunk_size,
                config.gcd_batch,
                config.segment_size,
            )
        finally:
            for item in reversed(values + owners):
                self.ffi.release(item)
        self._status(status)
        self._status(
            self.lib.b4_curve_options(
                self.context, int(self.atomic_actions), int(self.phase_clocks)
            )
        )
        if job["phase"] != "setup":
            self.plan = self._plan(job, context, budget)
            self.restore(job, budget)
        self.canonical_ready = True

    @staticmethod
    def _status(status):
        if status == -5:
            raise MemoryError("native curve preparation allocation failed")
        if status == -2:
            raise ValueError(
                "native curve table/metadata geometry is unsupported"
            )
        if status:
            raise ValueError("native canonical curve state refused")

    def feed(self, budget, *, _info=None):
        """Refill using a settled snapshot while the context is borrowed."""
        info = self._read() if _info is None else _info
        if not info.need:
            return
        if info.phase == 1 and self.plan is None:
            self.plan = self._plan(self.job, self.context_source, budget)
        lo = int(info.next)
        hi = min(int(info.hi), lo + 2 * self.config.segment_size)
        bound = self.job["b1"] if info.phase == 1 else None
        block, _ = self._block(lo, hi, bound, budget)
        self._status(self.lib.b4_curve_feed(self.context, block))

    def settle(self, budget):
        """Charge work and return a snapshot valid until context mutation."""
        info = self._read()
        work = int(info.work)
        delta = work - self.charged
        if delta < 0:
            raise ValueError("native curve work moved backwards")
        self.charged = work
        # This work has already happened, including discarded work. Record
        # every peer's debit before the root poll can stop publication.
        from ..execution.ecm_workers import _base_budget

        _base_budget(budget).record_incurred_work(delta)
        return info

    @exclusive_curve
    def advance_job(self, job, budget, context, config):
        poll_native_budget(budget, self)
        self.bind_job(job, budget, context, config)
        try:
            actions = (
                self.manifest.get("slice_actions", 128)
                if self.actions is None
                else self.actions
            )
            status = self.lib.b4_step(self.context, 5, actions)
            self.calls += 1
        finally:
            info = self.settle(budget)
        if status == -3:
            from ..execution.budget import BudgetExhaustedError

            raise BudgetExhaustedError("cancelled")
        if status < 0:
            self._status(status)
        if status == 2:
            self.feed(budget, _info=info)
            # Feed advances the cursor. Atomic export needs that new state.
            info = None
        poll_native_budget(budget, self)
        if status == 1:
            self.sync(
                allow_result=True,
                compact=not self.atomic_actions,
                _info=info,
            )
        elif self.atomic_actions:
            self.sync(allow_result=False, _info=info)

    def sync(self, *, allow_result=False, compact=False, _info=None):
        if (
            self.job is None
            or self.context is None
            or not self.canonical_ready
        ):
            return
        info = self._read() if _info is None else _info
        job = self.job
        first_result = bool(allow_result and info.done and not self.published)
        phase = PHASES[int(info.phase)]
        job["phase"] = phase
        if allow_result and info.done and compact:
            # Completed curves publish only checked identity/outcome metadata.
            # Tables stay private unless a diagnostic explicitly requests them.
            job["sigma"] = random.Random(job["seed"]).randrange(6, 2**63)
            factor = self._field(6)
            if factor and not utils.valid_divisor(factor, job["n"]):
                raise ValueError("native whole ECM returned an invalid factor")
            job["factor"], job["done"] = factor or None, True
            self.published = True
            if info.chain_chunks or info.chain and info.strict_state:
                job.update(
                    chain_identity=self.plan.identity,
                    chain_chunks=int(info.chain_chunks)
                    + int(bool(info.chain and info.strict_state)),
                    chain_replays=int(info.chain_replays)
                    + int(bool(info.chain and info.strict_state)),
                )
            for key in self.statistics:
                self.statistics[key] = int(getattr(info, key))
            self.phase_times = [int(v) for v in info.phase_ns]
            if first_result:
                for key, value in self.statistics.items():
                    self.lifetime_statistics[key] += value
                for index, value in enumerate(self.phase_times):
                    self.lifetime_phase_times[index] += value
            return
        cursor = dict(
            left=int(info.left),
            next=int(info.next),
            hi=int(info.hi),
            values=[],
            index=int(info.index),
        )
        bound = (
            job["b1"] if phase in ("setup", "stage_one", "replay") else None
        )
        entry = self.cursor_blocks.get((cursor["left"], cursor["next"], bound))
        if entry is not None:
            cursor["values"] = list(entry[1])
        if cursor["index"] > len(cursor["values"]):
            raise ValueError(
                "native checkpoint cursor lost its immutable block"
            )
        job["cursor"] = cursor
        powers = []
        for index in range(int(info.powers_count)):
            if self.lib.b4_curve_power(self.context, index, self.pair):
                raise ValueError("native checkpoint lost a prime power")
            powers.append([int(self.pair[0]), int(self.pair[1])])
        job["powers"] = powers
        job["product"] = self._field(5)
        job["terms"] = [
            self._field(9, index) for index in range(int(info.terms_count))
        ]
        job["recovery"] = int(info.recovery)
        if info.setup_done or (allow_result and info.done):
            job["sigma"] = random.Random(job["seed"]).randrange(6, 2**63)
        if info.setup_done:
            job["a24"] = self._field(0)
            job["value"] = [self._field(1), self._field(2)]
        if phase == "replay":
            job.update(
                replay_index=int(info.replay_index),
                replay_power=int(info.replay_power),
                replay_value=[self._field(3), self._field(4)],
            )
            if allow_result and info.done and info.operation == 7:
                job["replay_value"] = [self._field(14), self._field(15)]
                job["replay_power"] *= powers[job["replay_index"]][0]
        if phase in ("baby_steps", "giant_setup", "stage_two", "term_replay"):
            job["distance"] = int(info.distance)
            if info.distance >= 2:
                job["center"] = int(info.center)
                job["baby"] = [None] + [
                    [self._field(7, index), self._field(8, index)]
                    for index in range(1, int(info.baby_count) + 1)
                ]
                if phase in ("stage_two", "term_replay"):
                    job["giant"] = [self._field(10), self._field(11)]
                    job["previous"] = [self._field(12), self._field(13)]
        if info.chain_chunks:
            job.update(
                chain_identity=self.plan.identity,
                chain_chunks=int(info.chain_chunks),
                chain_replays=int(info.chain_replays),
            )
        if allow_result and info.done:
            factor = self._field(6)
            if factor and not utils.valid_divisor(factor, job["n"]):
                raise ValueError("native whole ECM returned an invalid factor")
            if info.chain and info.strict_state:
                job.update(
                    chain_identity=self.plan.identity,
                    chain_chunks=int(info.chain_chunks) + 1,
                    chain_replays=int(info.chain_replays) + 1,
                )
            if phase == "term_replay" and info.recovery < info.terms_count:
                job["recovery"] += 1
            job["factor"], job["done"] = factor or None, True
            self.published = True
        for key in self.statistics:
            self.statistics[key] = int(getattr(info, key))
        self.phase_times = [int(v) for v in info.phase_ns]
        if first_result:
            for key, value in self.statistics.items():
                self.lifetime_statistics[key] += value
            for index, value in enumerate(self.phase_times):
                self.lifetime_phase_times[index] += value

    def restore(self, job, budget):
        ffi = self.ffi
        info = ffi.new("b4_curve_info *")
        info.phase = PHASES.index(job["phase"])
        info.b1, info.b2 = job["b1"], job["b2"]
        info.setup_done = int(job["phase"] != "setup")
        info.chain_chunks, info.chain_replays = (
            job.get("chain_chunks", 0),
            job.get("chain_replays", 0),
        )
        for field, key in (
            ("center", "center"),
            ("distance", "distance"),
            ("replay_power", "replay_power"),
            ("replay_index", "replay_index"),
            ("recovery", "recovery"),
        ):
            setattr(info, field, job.get(key, 0))
        cursor = job["cursor"]
        info.left, info.next, info.hi, info.index = (
            cursor["left"],
            cursor["next"],
            cursor["hi"],
            cursor["index"],
        )
        info.powers_count, info.terms_count = (
            len(job["powers"]),
            len(job["terms"]),
        )
        info.baby_count = max(0, len(job.get("baby", [])) - 1)
        info.chain = int(self.plan is not None)
        bound = job["b1"] if job["phase"] in ("stage_one", "replay") else None
        block = ffi.NULL
        if cursor["next"] > cursor["left"]:
            block, _ = self._block(
                cursor["left"],
                cursor["next"],
                bound,
                budget,
                values=cursor["values"],
            )
        powers_block = ffi.NULL
        if job["powers"]:
            prime_values = [p for p, _ in job["powers"]]
            powers_block, _ = self._block(
                prime_values[0],
                prime_values[-1] + 1,
                job["b1"],
                budget,
                values=prime_values,
                purpose="powers",
            )
        values = [
            job.get("a24", 0),
            *job.get("value", [0, 0]),
            *job.get("replay_value", [0, 0]),
            job["product"],
            0,
            *job.get("giant", [0, 0]),
            *job.get("previous", [0, 0]),
        ]
        baby = [x for point in job.get("baby", [])[1:] for x in point]
        owners = []

        def pack(items):
            records = ffi.new("b4_value[]", max(1, len(items)))
            owners.append(records)
            for target, value in zip(records, items):
                raw = _bytes(value)
                data = ffi.new("unsigned char[]", raw)
                owners.append(data)
                target.data, target.size = data, len(raw)
            return records

        try:
            fields, points, terms = (
                pack(values),
                pack(baby),
                pack(job["terms"]),
            )
            self._status(
                self.lib.b4_curve_restore(
                    self.context,
                    info,
                    fields,
                    points,
                    terms,
                    powers_block,
                    block,
                )
            )
        finally:
            for owner in reversed(owners):
                ffi.release(owner)
            ffi.release(info)

    def discard(self, *, sync=True):
        if self.job is not None:
            try:
                if sync and not self.lib.b4_curve_info_get(
                    self.context, self.info
                ):
                    self.sync(allow_result=False)
            finally:
                self.lib.b4_abort(self.context)
                self.job = self.context_source = self.config = None
                self.canonical_ready = False
                self.cursor_blocks.clear()

    @exclusive_curve
    def close(self):
        try:
            if not self.closed:
                self.discard(sync=not self.published)
        finally:
            if not self.closed:
                close_resource_checks(self)
                if self.owns_blocks:
                    self.blocks.close()
                for owner in (self.info, self.atom, self.pair):
                    if owner is not None:
                        self.ffi.release(owner)
                self.info = self.atom = self.pair = None
            super().close()
