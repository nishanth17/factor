"""Exact coordinate permutations and translated adjacent Gray faces."""

from dataclasses import dataclass

from ..common import utils
from .families import MAX_A_FACTORS, _checksum

SCHEMA = "c10a07-sign-coordinates-v1"
ARMS = ("control", "complement_high", "swap_high")


@dataclass(frozen=True)
class SignSchedule:
    """Map Gray bits to sorted A factors; factor zero stays positive."""

    factor_count: int
    variable_order: tuple
    negative_mask: int = 0

    def __post_init__(self):
        utils.require_integer(self.factor_count, "factor_count", 1)
        if self.factor_count > MAX_A_FACTORS:
            raise ValueError("schedule exceeds the A-factor allowance")
        if not isinstance(self.variable_order, tuple) or (
            len(self.variable_order) != self.factor_count - 1
        ):
            raise ValueError("schedule order has the wrong shape")
        for index in self.variable_order:
            utils.require_integer(index, "factor coordinate", 1)
        if sorted(self.variable_order) != list(range(1, self.factor_count)):
            raise ValueError(
                "schedule order must permute variable coordinates"
            )
        utils.require_integer(self.negative_mask, "negative_mask", 0)
        if (
            self.negative_mask >= 1 << self.factor_count
            or self.negative_mask & 1
        ):
            raise ValueError(
                "initial signs must preserve fixed coordinate zero"
            )

    @property
    def canonical(self):
        return self.variable_order == tuple(range(1, self.factor_count)) and (
            self.negative_mask == 0
        )

    def mapping(self):
        return dict(
            schema=SCHEMA,
            factor_count=self.factor_count,
            fixed_factor_index=0,
            variable_order=list(self.variable_order),
            negative_mask=self.negative_mask,
        )

    @property
    def identity(self):
        return _checksum(self.mapping())

    @classmethod
    def from_mapping(cls, mapping):
        if not isinstance(mapping, dict) or set(mapping) != {
            "schema",
            "factor_count",
            "fixed_factor_index",
            "variable_order",
            "negative_mask",
        }:
            raise ValueError("invalid sign schedule mapping")
        if (
            mapping["schema"] != SCHEMA
            or type(mapping["fixed_factor_index"]) is not int
            or mapping["fixed_factor_index"] != 0
        ):
            raise ValueError("unknown sign schedule convention")
        order = mapping["variable_order"]
        if not isinstance(order, list) or len(order) > MAX_A_FACTORS - 1:
            raise ValueError("invalid sign schedule order")
        return cls(
            mapping["factor_count"], tuple(order), mapping["negative_mask"]
        )

    def mask(self, ordinal):
        """Return negative signs in physical, sorted-factor coordinates."""
        utils.require_integer(ordinal, "schedule ordinal", 0)
        if ordinal >= 1 << (self.factor_count - 1):
            raise ValueError("schedule ordinal is exhausted")
        gray = ordinal ^ (ordinal >> 1)
        mask = self.negative_mask
        for bit, coordinate in enumerate(self.variable_order):
            if gray & (1 << bit):
                mask ^= 1 << coordinate
        return mask

    def coverage(self, quota):
        utils.require_integer(quota, "quota", 1)
        if quota > min(8192, 1 << (self.factor_count - 1)):
            raise ValueError(
                "analytical coverage census exceeds its allowance"
            )
        masks = [self.mask(index) for index in range(quota)]
        return dict(
            schedule=self.mapping(),
            identity=self.identity,
            quota=quota,
            masks=masks,
            negative_counts=[
                sum(bool(m & (1 << i)) for m in masks)
                for i in range(self.factor_count)
            ],
        )


def schedule_for(factor_count, arm):
    """Return one of three declared arms without sweeping permutations."""
    order = tuple(range(1, factor_count))
    if arm == "control":
        return SignSchedule(factor_count, order)
    if arm == "complement_high" and factor_count >= 2:
        return SignSchedule(factor_count, order, 1 << (factor_count - 1))
    if arm == "swap_high" and factor_count >= 3:
        return SignSchedule(factor_count, order[:-2] + order[-2:][::-1])
    raise ValueError("unknown or inapplicable sign schedule")
