"""Bounded checked reuse, exact congruence indexing and useful-row feedback."""

from collections import OrderedDict

from ..common import prime_sieve, utils
from ..common.arithmetic import gcd, pow
from ..execution.resources import removed_arguments
from .factor_witness import factor_snapshot
from .indexed_relations import FactorRecord, encode_factors
from .lazy_relations import LazyCombinedRelation
from .relations import (
    AtomicRelation,
    CombinationResult,
    CombinedRelation,
    _combination_workspace,
    _combined_values,
    _verify_atomic,
    _verify_combined_values,
    parity_bits,
)


def payload_snapshot(row):
    """Compare complete immutable fields, including mutable-bypass faults.

    Packed bytes are copied into the witness; an ID, digest or a `verified`
    attribute is never sufficient to reuse a previous arithmetic check.
    """
    if isinstance(row, LazyCombinedRelation):
        return row.structural_snapshot()
    exponents = factor_snapshot(row.exponents)
    if isinstance(row, AtomicRelation):
        poly = row.polynomial
        return (
            row.relation_id,
            poly.n,
            poly.multiplier,
            poly.a,
            poly.b,
            poly.c,
            poly.square_coefficient,
            poly.identity,
            row.position,
            row.sign,
            exponents,
            row.residual,
            row.large_primes,
            getattr(row, "opaque", None),
        )
    return row.atom_ids, row.u, row.sign, exponents, row.square_correction


class CheckedStore:
    """Private proofs; restored or foreign objects start unchecked."""

    def __init__(self, base, config):
        self.base, self.config = base, config
        self.base_identity = (
            base.n,
            base.multiplier,
            base.bound,
            base.entries,
        )
        self.limit = config.checked_bytes
        self.records, self.primes, self.units = OrderedDict(), {}, {}
        self.used = 0
        self.canonical, self.canonical_ids = {}, {}
        self.canonical_used = 0
        self.canonical_limit = config.duplicate_bytes // 2

    def clear(self):
        self.records.clear()
        self.primes.clear()
        self.units.clear()
        self.canonical.clear()
        self.canonical_ids.clear()
        self.used = self.canonical_used = 0

    def _same_base(self):
        if self.base_identity != (
            self.base.n,
            self.base.multiplier,
            self.base.bound,
            self.base.entries,
        ):
            raise ValueError("checked store factor base changed")

    def prime(self, value, budget):
        self._same_base()
        reuse = self.config.recovery_policy == "checked-v2"
        if reuse and value in self.primes:
            budget.consume(1)
            return self.primes[value]
        budget.consume(value.bit_length() ** 2)
        result = utils.classify_prime(value) is utils.Primality.PROVEN
        # A small independent ceiling prevents proofs for one-shot rejects
        # from taking the whole row/atom attestation allowance.
        if reuse and len(self.primes) < 1024 and self.used + 512 <= self.limit:
            self.primes[value] = result
            self.used += 512
        return result

    def unit(self, value, budget):
        self._same_base()
        reuse = self.config.recovery_policy == "checked-v2"
        if reuse and value in self.units:
            budget.consume(1)
            return self.units[value]
        budget.consume(self.base.n.bit_length())
        result = gcd(value, self.base.n)
        if reuse and len(self.units) < 1024 and self.used + 512 <= self.limit:
            self.units[value] = result
            self.used += 512
        return result

    def _record(self, row, bits=None, atoms=()):
        if self.config.recovery_policy != "checked-v2":
            return
        reserve = 2048 + 128 * len(row.exponents)
        if isinstance(row, CombinedRelation):
            reserve += 128 * len(row.atom_ids)
            reserve += sum(2048 + 128 * len(a.exponents) for a in atoms)
        if bits is not None:
            reserve += 8 * ((bits.bit_length() + 7) // 8)
        old = self.records.pop(id(row), None)
        if old is not None:
            self.used -= old[3]
        # Keep recent published work useful after the store outgrows this
        # finite cache. An evicted proof never evicts its retained relation.
        while self.records and self.used + reserve > self.limit:
            _, retired = self.records.popitem(last=False)
            self.used -= retired[3]
        if self.used + reserve <= self.limit:
            snapshot = payload_snapshot(row)
            witnesses = tuple((a, payload_snapshot(a)) for a in atoms)
            self.records[id(row)] = row, snapshot, bits, reserve, witnesses
            self.used += reserve

    def _known(self, row, budget):
        self._same_base()
        record = self.records.get(id(row))
        budget.consume(len(row.exponents) + 1)
        if record is not None and record[0] is row:
            if record[1] != payload_snapshot(row):
                raise ValueError("checked immutable relation was modified")
            self.records.move_to_end(id(row))
            return record
        return None

    def atom(self, atom, budget):
        if self.config.recovery_policy == "checked-v2" and self._known(
            atom, budget
        ):
            return
        _verify_atomic(
            atom,
            self.base,
            residual_bound=self.config.residual_bound,
            large_prime_bound=self.config.large_prime_bound,
            large_product_bound=self.config.large_product_bound,
            small_prime_bound=self.config.small_prime_bound,
            residual_arity=getattr(self.config, "residual_arity", 2),
            opaque_pairs=getattr(self.config, "opaque_pairs", False),
            budget=budget,
            prime_check=(
                self.prime
                if self.config.recovery_policy == "checked-v2"
                else None
            ),
        )
        self._record(atom)

    @removed_arguments(("memory_bytes",))
    def combine(self, atoms, budget):
        identities = tuple(atom.relation_id for atom in atoms)
        if not 1 <= len(atoms) <= 256 or len(set(identities)) != len(atoms):
            raise ValueError("invalid checked combination provenance")
        _combination_workspace(atoms, self.base)
        for atom in atoms:
            self.atom(atom, budget)
            for prime in atom.residual_primes:
                divisor = self.unit(prime, budget)
                if utils.valid_divisor(divisor, self.base.n):
                    return CombinationResult(None, divisor)
                if divisor != 1:
                    raise ValueError("nonunit residual has no proper split")
        if getattr(self.config, "lazy_cycles", False):
            row = LazyCombinedRelation(atoms, self.base.n, budget)
        else:
            budget.consume(
                sum(len(a.exponents) + 1 for a in atoms)
                * (self.base.n.bit_length() + 1)
            )
            row = CombinedRelation(
                identities, *_combined_values(atoms, self.base.n)
            )
        # Exact atomic identities and even residual multiplicities establish
        # the combined identity. Public verification stays independent.
        self._record(row, atoms=atoms)
        return CombinationResult(row, None)

    @removed_arguments(("memory_bytes",))
    def row(self, row, store, budget):
        if (
            isinstance(row, LazyCombinedRelation)
            and row._modulus != self.base.n
        ):
            raise ValueError("lazy relation belongs to a different target")
        record = (
            self._known(row, budget)
            if self.config.recovery_policy == "checked-v2"
            else None
        )
        if isinstance(row, AtomicRelation):
            if row.residual_product != 1:
                raise ValueError("partial atom is not a matrix row")
            if record is None:
                self.atom(row, budget)
        else:
            atoms = tuple(store.get(i) for i in row.atom_ids)
            budget.consume(
                sum(len(a.exponents) + 1 for a in atoms if a is not None)
            )
            if record is not None and (
                len(record[4]) != len(atoms)
                or any(
                    a is not witness[0] or payload_snapshot(a) != witness[1]
                    for a, witness in zip(atoms, record[4])
                )
            ):
                record = None
            if record is None:
                for identity in row.atom_ids:
                    atom = store.get(identity)
                    if atom is None or atom.relation_id != identity:
                        raise ValueError("missing checked provenance")
                    self.atom(atom, budget)
                _combination_workspace(atoms, self.base)
                if isinstance(row, LazyCombinedRelation):
                    row.verify_structure(atoms, budget)
                else:
                    _verify_combined_values(row, self.base, atoms, budget)
            if isinstance(row, LazyCombinedRelation):
                row.bind_budget(budget)
        if record is not None and record[2] is not None:
            return record[2]
        bits = parity_bits(row, self.base)
        self._record(
            row, bits, () if isinstance(row, AtomicRelation) else atoms
        )
        return bits

    def congruence(self, atom, budget):
        """Return an exact key or a divisor before any duplicate is dropped."""
        correction = atom.square_correction
        if correction == 1:
            self._same_base()
            divisor = 1
        else:
            divisor = self.unit(correction, budget)
        if utils.valid_divisor(divisor, self.base.n):
            return None, divisor
        if divisor != 1:
            return None, None
        budget.consume(self.base.n.bit_length() ** 2 + len(atom.exponents))
        root = (
            atom.u % self.base.n
            if correction == 1
            else atom.u * pow(correction, -1, self.base.n) % self.base.n
        )
        factors = atom.exponents
        encoded = (
            factors.encoded()
            if type(factors) is FactorRecord
            and factors._primes is self.base.primes
            else encode_factors(factors, self.base)
        )
        key = (
            min(root, self.base.n - root),
            atom.sign,
            encoded,
            atom.residual_primes,
        )
        return key, None

    def duplicate(self, key, store, budget):
        old = self.canonical.get(key)
        if old is None or store.get(old.relation_id) is not old:
            return False
        # The canonical key is exact, but its representative may have been
        # corrupted through an unsupported mutation bypass; never trust it.
        self.atom(old, budget)
        checked, divisor = self.congruence(old, budget)
        return divisor is None and checked == key

    def remember(self, atom, row=None, key=None, store=None):
        """Publish only inside the collector's allocation-failure boundary."""
        self._record(atom)
        if row is not None:
            self._record(row, atoms=tuple(store[i] for i in row.atom_ids))
        if key is not None and key not in self.canonical:
            reserve = 2048 + 256 * len(atom.exponents)
            if self.canonical_used + reserve <= self.canonical_limit:
                self.canonical[key] = atom
                self.canonical_ids[atom.relation_id] = key, reserve
                self.canonical_used += reserve

    def forget(self, atom):
        self.retire_proof(atom)
        entry = self.canonical_ids.pop(atom.relation_id, None)
        if entry is not None:
            key, reserve = entry
            self.canonical.pop(key, None)
            self.canonical_used -= reserve

    def retire_proof(self, row):
        """Release an obsolete witness without removing its successor."""
        record = self.records.pop(id(row), None)
        if record is not None:
            self.used -= record[3]


@removed_arguments(("memory_bytes",))
def residual_floor(base, budget):
    """First omitted eligible prime, conservatively, even for custom bases.

    Every prime dividing the remaining norm must divide U²-kN. Thus a prime
    below B can remain only if kN is a square modulo it and it was omitted.
    A custom incomplete base or multiplier exception lowers this bound.
    """
    budget.consume(base.bound)
    for prime in prime_sieve.prime_sieve(base.bound):
        if prime in base._columns:
            continue
        budget.consume(prime.bit_length() ** 2)
        residue = base.n_prime % prime
        if (
            residue == 0
            or prime == 2
            or pow(residue, (prime - 1) // 2, prime) == 1
        ):
            return prime
    return base.bound


class UsefulTarget:
    """Derived occupied columns and bounded filtered-deficit feedback."""

    def __init__(self):
        self.count = self.columns = self.next_target = 0
        self.last_raw = 0
        self.last_useful = 0

    def observe(self, rows, collector, budget):
        if self.count > len(rows):
            raise ValueError("useful-target cursor exceeds retained rows")
        while self.count < len(rows):
            row = rows[self.count]
            bits = collector._checked.row(row, collector._atoms, budget)
            budget.consume(1 + (bits.bit_length() + 63) // 64)
            self.columns |= bits
            self.count += 1

    def target(self, surplus):
        return self.next_target or self.columns.bit_count() + surplus

    def filtered(self, count, stats, surplus, growth, limit):
        useful = (
            stats["output_rows"]
            - stats["output_columns"]
            + stats["zero_dependencies"]
        )
        deficit = max(1, surplus - useful)
        incoming = max(1, count - self.last_raw)
        survived = max(1, useful - self.last_useful)
        desired = (deficit * incoming + survived - 1) // survived
        # Extrapolation is uncertain while the graph matures. A quarter of
        # the current row count bounds overshoot; the normal growth floor
        # prevents repeatedly preparing nearly identical inputs.
        topup = min(max(growth, (count + 3) // 4), max(growth, desired))
        self.next_target = (
            count + topup if limit is None else min(limit, count + topup)
        )
        self.last_raw, self.last_useful = count, useful

    def state(self):
        return dict(
            count=self.count,
            columns=hex(self.columns),
            next_target=self.next_target,
            last_raw=self.last_raw,
            last_useful=self.last_useful,
        )

    def restore(self, state, rows, collector, budget):
        if not isinstance(state, dict) or set(state) != set(self.state()):
            raise ValueError("invalid useful-target state")
        count = state["count"]
        if type(count) is not int or not 0 <= count <= len(rows):
            raise ValueError("invalid useful-target cursor")
        self.observe(rows[:count], collector, budget)
        if state["columns"] != hex(self.columns):
            raise ValueError("useful-target columns differ from checked rows")
        for name in ("next_target", "last_raw", "last_useful"):
            value = state[name]
            if type(value) is not int or not -len(
                collector.factor_base.entries
            ) - 1 <= value <= max(len(rows), state["next_target"]):
                raise ValueError("useful-target feedback exceeds its bound")
        if (
            not 0 <= state["last_raw"] <= count
            or not 0 <= state["next_target"]
            or state["last_useful"] > state["last_raw"]
        ):
            raise ValueError("invalid useful-target feedback order")
        self.next_target = state["next_target"]
        self.last_raw, self.last_useful = (
            state["last_raw"],
            state["last_useful"],
        )
