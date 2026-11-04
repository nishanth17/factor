"""Bounded single-process ECM parameter replay on explicit common curve seeds.

This diagnostic contract is separate from production Manifest seed derivation.
Every replay ID binds both the original curve and selected configuration. It
prices complete curves and setup; it makes no portfolio performance claim.
"""

import argparse
import fcntl
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path

REPLAY_VERSION = "b12-f4-common-curve-replay-v1"


def replay_assignment(reference, selected_identity):
    """Bind changed parameters without silently changing the curve seed."""
    from ....execution.worker_contract import digest

    result = dict(reference)
    result["id"] = digest([REPLAY_VERSION, reference["id"], selected_identity])
    return result


def child(path):
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    fd = int(os.environ["FACTOR_DIAGNOSTIC_LOCK_FD"])
    actual, expected = (
        os.fstat(fd),
        os.stat("/private/tmp/factor-performance.lock"),
    )
    if (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino):
        raise RuntimeError("lock identity changed")
    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    request = json.loads(path.read_text())
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_CPU, (40, 41))
    resource.setrlimit(resource.RLIMIT_FSIZE, (4 * 2**20, 4 * 2**20))
    from ....execution.ecm_worker import CurveEngine
    from ....execution.ecm_workers import arithmetic_options
    from ....execution.portfolio_config import PortfolioConfig
    from ....execution.worker_contract import Manifest, WorkerSettings

    settings = WorkerSettings(**request["settings"])
    base = PortfolioConfig(**request["config"])
    selected = PortfolioConfig(
        **dict(request["config"], **request["parameters"])
    )
    original = Manifest(
        request["n"],
        request["seed"],
        "root",
        arithmetic_options(base, settings),
    )
    plan = Manifest(
        request["n"],
        request["seed"],
        "root",
        arithmetic_options(selected, settings),
    )
    if original.tiers != plan.tiers:
        raise ValueError(
            "parameter replay must preserve numerical search bounds"
        )
    count = request["prefix"]
    if type(count) is not int or not 1 <= count <= min(16, len(original.jobs)):
        raise ValueError("finite prefix outside original manifest")
    engine = CurveEngine(settings)
    reference, calls, work = None, [], 0

    def invoke():
        nonlocal reference, work
        if time.monotonic() - request["entry"] >= request["wall"] - 2:
            raise RuntimeError("replay wall allowance exhausted")
        if time.process_time() >= request["cpu"] - 1:
            raise RuntimeError("replay CPU allowance exhausted")
        before = time.monotonic()
        rows = []
        for position in range(count):
            curve = original.assignment(position)
            assignment = replay_assignment(curve, plan.identity)
            task = dict(
                assignment=assignment,
                options=plan.options,
                schedule=plan.metadata,
                token=position + 1,
                attempt=1,
                stage="curve",
                work=min(32_000_000, 10**12 - work),
                cpu=min(2, request["cpu"] - time.process_time()),
                deadline=request["entry"] + request["wall"] - 2,
            )
            result = engine.execute(task, cancelled=lambda: False)
            work += result["work"]
            if (
                result["reason"] != "complete"
                or result["sigma"] != curve["sigma"]
            ):
                raise ValueError("incomplete or changed replay curve")
            divisor = result["factor"]
            if divisor is not None and not (
                1 < divisor < request["n"] and request["n"] % divisor == 0
            ):
                raise ValueError("nonproper replay divisor")
            rows.append(
                [
                    curve["id"],
                    result["id"],
                    result["sigma"],
                    divisor,
                    result["phase"],
                    result["chain_chunks"],
                ]
            )
        if reference is None:
            reference = rows
        elif reference != rows:
            raise ValueError("replay result changed between repetitions")
        return time.monotonic() - before

    started = time.monotonic()
    warmup_calls = 0
    for _ in range(2048):
        invoke()
        warmup_calls += 1
        if time.monotonic() - started >= 3:
            break
    else:
        raise RuntimeError("validated warmup did not reach three seconds")
    warmup = time.monotonic() - started
    for _ in range(request["repetitions"]):
        calls.append(invoke())
    result = dict(
        complete=True,
        contract=REPLAY_VERSION,
        rows=reference,
        calls=calls,
        warmup_seconds=warmup,
        warmup_calls=warmup_calls,
        work=work,
        cpu=time.process_time(),
        preparations=engine.preparations,
        decode_hits=engine.decode_hits,
        decode_misses=engine.decode_misses,
        decode_refusals=engine.decode_refusals,
        original_manifest=original.identity,
        selected_manifest=plan.identity,
    )
    capture = Path(request["capture"])
    capture.mkdir(exist_ok=False)
    with (capture / "result.json").open("x") as stream:
        json.dump(result, stream, allow_nan=False)
    print(
        json.dumps(
            dict(
                complete=True,
                calls=calls,
                signature=hashlib.sha256(
                    json.dumps(reference).encode()
                ).hexdigest(),
            )
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    child(parser.parse_args().request)


if __name__ == "__main__":
    main()
