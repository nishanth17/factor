"""Build a versioned selector catalogue from reviewed literal costs."""

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from ....execution.ecm_selection import (
    COMPATIBILITY,
    PROFILE_BYTES,
    PROFILE_VERSION,
    SelectionCosts,
)
from ....execution.ecm_workers import arithmetic_options
from ....execution.portfolio_config import PortfolioConfig
from ....execution.worker_contract import (
    SEARCH_VERSION,
    VERSION,
    WorkerSettings,
    canonical,
    digest,
    integer,
)


def build_profile(reviewed, *, status="training"):
    """Record costs without fitting them or inspecting any factor labels."""
    if status not in ("training", "frozen-for-confirmation", "confirmed"):
        raise ValueError("unsupported catalogue status")
    if set(reviewed) != {"serial_settings", "rows"}:
        raise ValueError("reviewed input must contain settings and rows")
    if not isinstance(reviewed["rows"], list) or len(reviewed["rows"]) > 256:
        raise ValueError("finite reviewed row set required")
    rows, keys = [], set()
    serial = WorkerSettings(**reviewed["serial_settings"])
    if serial.workers != 1 or serial.worker_memory != 64 * 2**20:
        raise ValueError("incompatible serial mathematical context")
    for row in reviewed["rows"]:
        if set(row) != {"digits", "config", "cold", "warm"}:
            raise ValueError("unsupported reviewed fields")
        digits = integer(row["digits"], "digits", 2, 1234)
        config = PortfolioConfig(**row["config"])
        options = digest(arithmetic_options(config, serial))
        key = digits, options
        if key in keys:
            raise ValueError("duplicate reviewed regime")
        keys.add(key)
        selected = dict(digits=digits, options=options, tiers=config.ecm_tiers)
        for condition in ("cold", "warm"):
            candidate = row[condition]
            if set(candidate) != {"settings", "costs"}:
                raise ValueError("invalid reviewed candidate")
            settings = WorkerSettings(**candidate["settings"])
            if settings.worker_memory != serial.worker_memory:
                raise ValueError("mathematical memory differs between arms")
            if digest(arithmetic_options(config, settings)) != options:
                raise ValueError("candidate changes mathematical search")
            costs = SelectionCosts(
                digits=digits, options=options, **candidate["costs"]
            )
            selected[condition] = dict(
                settings=asdict(settings), costs=asdict(costs)
            )
        rows.append(selected)
    profile = dict(
        version=PROFILE_VERSION,
        status=status,
        executor_version=VERSION,
        search_version=SEARCH_VERSION,
        compatibility=COMPATIBILITY,
        serial_settings=asdict(serial),
        rows=rows,
    )
    if len(canonical(profile).encode("ascii")) + 1 > PROFILE_BYTES:
        raise ValueError("reviewed catalogue exceeds encoded ownership")
    return profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reviewed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--status",
        choices=("training", "frozen-for-confirmation", "confirmed"),
        default="training",
    )
    args = parser.parse_args()
    if args.reviewed.stat().st_size > PROFILE_BYTES:
        raise ValueError("reviewed input exceeds cap")
    profile = build_profile(
        json.loads(args.reviewed.read_text()), status=args.status
    )
    with args.output.open("x") as stream:
        stream.write(canonical(profile) + "\n")
    print(hashlib.sha256(args.output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
