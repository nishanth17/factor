"""Run the native acceptance suite and save every test's outcome."""

import argparse
import json
import time
import unittest
from pathlib import Path

from ..support.paths import (
    PACKAGE_ROOT,
    REPOSITORY_ROOT,
)
from .phase_one import environment


class RecordedResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.records = []

    def startTest(self, test):
        self.started = time.perf_counter()
        super().startTest(test)

    def _record(self, test, status):
        self.records.append(
            {
                "test": test.id(),
                "status": status,
                "seconds": time.perf_counter() - self.started,
            }
        )

    def addSuccess(self, test):
        self._record(test, "passed")
        super().addSuccess(test)

    def addFailure(self, test, err):
        self._record(test, "failed")
        super().addFailure(test, err)

    def addError(self, test, err):
        self._record(test, "error")
        super().addError(test, err)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--milestone", default="phase_1")
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.discover(
        str(PACKAGE_ROOT / "tests"),
        top_level_dir=str(REPOSITORY_ROOT),
    )
    runner = unittest.TextTestRunner(verbosity=2, resultclass=RecordedResult)
    start = time.perf_counter()

    result = runner.run(suite)
    args.output.write_text(
        json.dumps(
            {
                "milestone": args.milestone,
                "environment": environment(),
                "seconds": time.perf_counter() - start,
                "tests_run": result.testsRun,
                "successful": result.wasSuccessful(),
                "failures": len(result.failures),
                "errors": len(result.errors),
                "tests": result.records,
            },
            indent=2,
        )
        + "\n"
    )
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
