"""Read-only A1b terminal audit after the original analysis adapter failure.

This separate identity corrects generator materialization only. It preserves
all original timing, source, reservation and analysis failures and never runs
new factoring observations or changes the frozen advancement interpretation.
"""

import json
import time
from unittest.mock import patch

from . import a1b
from .c10a04.accounting import Ledger, digest, output_lock, receipt


def geometric_mean(values):
    """Consume a finite iterator once before applying the frozen formula."""
    values = tuple(values)
    return a1b.exp(sum(map(a1b.log, values)) / len(values))


def audit():
    """Reconcile receipts and apply the original formula without new spend."""
    started = time.perf_counter()
    protocol = a1b.load()
    with output_lock(a1b.OUTPUT):
        ledger = Ledger(a1b.OUTPUT, digest(a1b.FREEZE), a1b.limits(protocol))
        # This is a reporting/QA repair, not a new measurement or phase retry.
        # The failed 60-second analysis reservation stays charged in full.
        with patch.object(a1b, "geometric_mean", geometric_mean):
            result = a1b.analyze(protocol, ledger)
        calls = []
        for entry in ledger.data["calls"]:
            path = a1b.OUTPUT / entry["receipt"]
            if entry["status"] == "complete":
                if digest(path) != entry["sha256"]:
                    raise ValueError("terminal receipt mismatch")
                calls.append(json.loads(path.read_text()))
        result["identity"] = "a1b-terminal-audit-v1"
        result["repair"] = (
            "Original analysis failed on len(generator); materialize the "
            "finite input ratios without changing any analysis rule. "
            "Read-only audit under ordinary runner-regression/reporting QA; "
            "failed analysis reservation remains charged, no new calls."
        )
        result["source_freeze_sha256"] = digest(a1b.FREEZE)
        result["charges"] = {
            phase: {
                key: ledger.used(phase, key) for key in ("wall", "cpu", "work")
            }
            for phase in protocol["accounts"]
            if phase != "total"
        }
        result["total_charges"] = {
            key: ledger.used(None, key) for key in ("wall", "cpu", "work")
        }
        result["failures"] = [
            row for row in ledger.data["calls"] if row["status"] != "complete"
        ]
        actual = {key: 0.0 for key in ("wall", "cpu")}
        maxima = {key: 0 for key in ("rss", "owned", "queue")}
        for call in calls:
            row = call["result"]
            if "physical" in row:
                for key, name in (
                    ("wall", "actual_wall"),
                    ("cpu", "actual_cpu"),
                ):
                    actual[key] += call[name]
                maxima["rss"] = max(maxima["rss"], row["rss_bytes"])
                maxima["owned"] = max(
                    maxima["owned"], row["stats"].get("workspace_bytes", 0)
                )
                maxima["queue"] = max(maxima["queue"], row["queue_peak"])
        # The diagnostic's own durable import explicitly contains its clocks.
        diagnostic = next(
            c["result"]
            for c in calls
            if c["identity"].get("kind") == "imported diagnostic"
        )
        for key in actual:
            actual[key] += diagnostic["actual_" + key]
        result["factoring_and_diagnostic_actual"] = actual
        result["maxima_bytes"] = maxima
        result["receipt_count"] = len(calls)
        result["audit_seconds"] = time.perf_counter() - started
        return result


def main():
    result = audit()
    receipt(a1b.OUTPUT / "terminal-audit-v1.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
