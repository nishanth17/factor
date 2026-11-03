"""Queued pure scheduler QA under the released shared performance window.

This launcher cannot authorize execution while C9a owns timing. Preparation
phase used neither deterministic-harness opportunity to execute this code.
"""

import argparse
import fcntl
import hashlib
import json
import math
import resource
import signal
import sys
import time
import unittest
from pathlib import Path


def atomic(path, value):
    temporary = path.with_suffix(".pending")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--window-released", action="store_true")
    parser.add_argument("--receipts", type=Path, required=True)
    args = parser.parse_args()
    if not args.window_released:
        parser.error("requires explicit release of C9a's performance window")
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    root = Path(__file__).resolve().parents[4]
    pins = {
        str(path.relative_to(root)): hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in (
            root / "v2/execution/worker_contract.py",
            root / "v2/tests/test_worker_contract.py",
            Path(__file__).resolve(),
        )
    }
    args.receipts.mkdir(parents=True, exist_ok=True)
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = args.receipts / "pure-ledger.json"
        ledger = (
            json.loads(path.read_text())
            if path.exists()
            else dict(wall=0, cpu=0, attempts=[])
        )
        remaining = min(60 - ledger["wall"], 60 - ledger["cpu"])
        if remaining < 1:
            raise RuntimeError("cumulative pure-harness allowance exhausted")
        receipt = args.receipts / f"pure-{len(ledger['attempts']) + 1}.json"
        # Full reservation survives interruption; only this completed receipt
        # can demonstrate unused resources and refund them.
        ledger["wall"] += remaining
        ledger["cpu"] += remaining
        ledger["attempts"].append(str(receipt))
        atomic(path, ledger)
        cpu_start, started = time.process_time(), time.monotonic()
        cpu_cap = max(1, math.floor(cpu_start + remaining))
        resource.setrlimit(resource.RLIMIT_CPU, (cpu_cap, cpu_cap))

        def timeout(*_):
            raise KeyboardInterrupt("pure harness wall allowance")

        signal.signal(signal.SIGALRM, timeout)
        signal.setitimer(signal.ITIMER_REAL, remaining)
        result, error = None, None
        try:
            suite = unittest.defaultTestLoader.loadTestsFromName(
                "v2.tests.test_worker_contract"
            )
            result = unittest.TextTestRunner(verbosity=2).run(suite)
            forbidden = (
                "v2.ecm",
                "v2.portfolio",
                "v2.execution.ecm_worker",
                "v2.execution.ecm_workers",
            )
            if any(
                name == prefix or name.startswith(prefix + ".")
                for name in sys.modules
                for prefix in forbidden
            ):
                raise RuntimeError(
                    "pure harness imported arithmetic/worker runtime"
                )
        except BaseException as failure:
            error = repr(failure)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        wall, cpu = time.monotonic() - started, time.process_time() - cpu_start
        payload = dict(
            source_sha256=pins,
            wall=wall,
            cpu=cpu,
            passed=bool(result and result.wasSuccessful() and error is None),
            tests=None if result is None else result.testsRun,
            error=error,
            children_started=0,
            arithmetic_executed=False,
        )
        with receipt.open("x") as output:
            json.dump(payload, output, indent=2, allow_nan=False)
        ledger["wall"] -= max(0, remaining - wall)
        ledger["cpu"] -= max(0, remaining - cpu)
        atomic(path, ledger)
        if not payload["passed"]:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
