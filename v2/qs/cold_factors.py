"""Cheap-factor recomputation and bounded persistent cold records."""

import hashlib
import os
import tempfile
from collections import OrderedDict
from collections.abc import Sequence
from dataclasses import FrozenInstanceError

from .indexed_relations import decode_factors


class RecomputedFactors(Sequence):
    """Omit identified cheap factors, retaining exact provenance."""

    def __init__(self, polynomial, position, small, retained, count):
        for name, value in (
            ("polynomial", polynomial),
            ("position", position),
            ("small", tuple(small)),
            ("retained", retained),
            ("count", count),
        ):
            object.__setattr__(self, name, value)

    def __setattr__(self, name, value):
        raise FrozenInstanceError("recomputed factors are immutable")

    def __delattr__(self, name):
        raise FrozenInstanceError("recomputed factors are immutable")

    def __len__(self):
        return self.count

    def __iter__(self):
        value = abs(self.polynomial.a * self.polynomial.value(self.position))
        correction = self.polynomial.square_coefficient
        if correction == 1:
            remainder = 0
        else:
            value, remainder = divmod(value, correction**2)
        if not value or remainder:
            raise ValueError("invalid recomputation provenance")
        for prime in self.small:
            exponent = 0
            while value % prime == 0:
                value //= prime
                exponent += 1
            if not exponent:
                raise ValueError("omitted factor no longer divides its norm")
            yield prime, exponent
        yield from self.retained

    def __getitem__(self, index):
        return tuple(self)[index]

    def __eq__(self, other):
        return tuple(self) == tuple(other)

    def __hash__(self):
        return hash(tuple(self))


class SpillStore:
    """Append-only owned file with bounded LRU reads and explicit cleanup.

    Files persist for checkpoints. No destructor deletes evidence. References
    bind filename, extent, checksum, factor base and canonical record shape.
    A refused append never publishes an offset; partial files remain owned.
    """

    def __init__(self, directory, limit, cache_bytes, budget):
        self.directory = os.path.realpath(directory)
        self.limit, self.cache_limit, self.budget = limit, cache_bytes, budget
        self.used, self.cache_used = 0, 0
        self.cache = OrderedDict()
        self.files = {}
        self.writer = None

    def _path(self, name):
        if (
            not isinstance(name, str)
            or not name.startswith("c9a-cold-")
            or (os.path.basename(name) != name or not name.endswith(".bin"))
        ):
            raise ValueError("invalid owned spill filename")
        return os.path.join(self.directory, name)

    def append(self, data):
        if len(data) + self.used > self.limit:
            raise MemoryError("production disk allowance exhausted")
        self.budget.consume(len(data) + 1)
        if self.writer is None:
            if len(self.files) >= 1024:
                raise MemoryError(
                    "production spill file-count allowance exhausted"
                )
            fd, path = tempfile.mkstemp(
                prefix="c9a-cold-", suffix=".bin", dir=self.directory
            )
            self.writer = os.fdopen(fd, "ab", buffering=0)
            self.name = os.path.basename(path)
            self.files[self.name] = 0
        offset = self.files[self.name]
        try:
            written = self.writer.write(data)
        except OSError:
            extent = os.fstat(self.writer.fileno()).st_size
            self.used += max(0, extent - offset)
            self.files[self.name] = extent
            raise
        if written != len(data):
            self.used += max(0, written or 0)
            self.files[self.name] += max(0, written or 0)
            raise OSError("partial cold-record append")
        self.files[self.name] += len(data)
        self.used += len(data)
        return self.name, offset, len(data), hashlib.sha256(data).hexdigest()

    def read(self, reference):
        if not isinstance(reference, (tuple, list)) or len(reference) != 4:
            raise ValueError("invalid cold-record reference")
        name, offset, length, digest = reference
        if (
            not isinstance(name, str)
            or not isinstance(digest, str)
            or len(digest) != 64
            or type(offset) is not int
            or type(length) is not int
            or (
                offset < 0
                or not 0 <= length <= self.cache_limit
                or offset + length > self.files.get(name, -1)
            )
        ):
            raise ValueError("spill extent exceeds its owned bound")
        self.budget.consume(length + 1)
        # Always check file bytes: a stale cache must not hide external damage.
        fd = os.open(
            self._path(name), os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
        )
        with os.fdopen(fd, "rb") as stream:
            stream.seek(offset)
            data = stream.read(length)
        if len(data) != length or hashlib.sha256(data).hexdigest() != digest:
            raise ValueError("missing, truncated or corrupt cold record")
        return data

    def decoded(self, reference, primes, count):
        """Reuse bounded decoding only after checking current file bytes."""
        data, key = self.read(reference), tuple(reference)
        previous = self.cache.pop(key, None)
        if previous is not None:
            self.cache_used -= previous[4]
            if (
                previous[0] == data
                and previous[1] is primes
                and previous[2] == count
            ):
                self.cache[key] = previous
                self.cache_used += previous[4]
                return iter(previous[3])
        reserve = 128 + len(data) + 192 * count
        if reserve > self.cache_limit:
            return decode_factors(data, primes, count)
        while self.cache and self.cache_used + reserve > self.cache_limit:
            _, retired = self.cache.popitem(last=False)
            self.cache_used -= retired[4]
        self.budget.consume(count + 1)
        pairs = tuple(decode_factors(data, primes, count))
        self.cache[key] = data, primes, count, pairs, reserve
        self.cache_used += reserve
        return iter(pairs)

    def manifest(self):
        if self.writer is not None:
            self.writer.flush()
            os.fsync(self.writer.fileno())
        return [[name, size] for name, size in sorted(self.files.items())]

    def restore_manifest(self, manifest):
        if (
            self.files
            or not isinstance(manifest, list)
            or len(manifest) > 1024
        ):
            raise ValueError("invalid spill manifest")
        files, total = {}, 0
        for name, size in manifest:
            self._path(name)
            if name in files or type(size) is not int or size < 0:
                raise ValueError("invalid spill extent")
            total += size
            if total > self.limit:
                raise MemoryError("saved production disk allowance exhausted")
            if (
                os.path.islink(self._path(name))
                or os.stat(self._path(name)).st_size < size
            ):
                raise ValueError("missing or truncated spill file")
            files[name] = size
        self.files, self.used = files, total

    def close(self):
        if self.writer is not None:
            self.writer.close()
            self.writer = None

    def cleanup(self):
        """Explicitly invalidate this store and its external checkpoints."""
        self.close()
        for name in self.files:
            os.unlink(self._path(name))
        self.files.clear()
        self.cache.clear()
        self.used = self.cache_used = 0


class ColdFactors(Sequence):
    def __init__(self, owner, reference, primes, count):
        if (
            type(count) is not int
            or not 0 <= count <= len(primes)
            or not isinstance(reference, (tuple, list))
            or len(reference) != 4
            or type(reference[2]) is not int
            or not 0 <= reference[2] <= 10 * count
        ):
            raise ValueError("cold encoding exceeds its declared factor count")
        for name, value in (
            ("owner", owner),
            ("reference", tuple(reference)),
            ("primes", primes),
            ("count", count),
        ):
            object.__setattr__(self, name, value)

    def __setattr__(self, name, value):
        raise FrozenInstanceError("cold factor references are immutable")

    def __delattr__(self, name):
        raise FrozenInstanceError("cold factor references are immutable")

    def encoded(self):
        return self.owner.read(self.reference)

    def __len__(self):
        return self.count

    def __iter__(self):
        return self.owner.decoded(self.reference, self.primes, self.count)

    def __getitem__(self, index):
        return tuple(self)[index]

    def __eq__(self, other):
        return tuple(self) == tuple(other)

    def __hash__(self):
        return hash(tuple(self))
