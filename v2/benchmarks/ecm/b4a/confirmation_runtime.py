"""Finite ECM confirmation protocol with complete retained-owner accounting."""

import json
import sys
import time


def close_spawn_helper():
    """Reap this isolated interpreter's helper before its parent waits."""
    module = sys.modules.get("multiprocessing.resource_tracker")
    if module is not None:
        tracker = module._resource_tracker
        if tracker._pid is not None:
            tracker._stop()


def measured_call(session, request):
    """Timestamp the existing checked publication boundary once per split."""
    from .... import portfolio

    original = portfolio._split
    publications = []

    def publish(state, divisor):
        original(state, divisor)
        publications.append(dict(at=time.monotonic(), divisor=int(divisor)))

    portfolio._split = publish
    started = time.monotonic()
    try:
        row = session.call(request)
    finally:
        portfolio._split = original
    row.update(
        dispatch_started=started,
        dispatch_wall=time.monotonic() - started,
        publications=publications,
    )
    return row


def serve(*, cold=False):
    from .production_batch import ProductionSession

    session = ProductionSession()
    warmup_failed = None
    try:
        for line in sys.stdin:
            request = json.loads(line)
            op = request.get("op", "call")
            if op == "warm":
                useful = 0.0
                calls = 0
                started = time.monotonic()
                while useful < 3 and calls < 5000:
                    latest = measured_call(session, request)
                    useful += latest["api_wall"]
                    calls += 1
                    if latest["resource_stop"] or latest.get(
                        "native_assertion_failure"
                    ):
                        warmup_failed = latest
                        break
                    if time.monotonic() - started > 45:
                        warmup_failed = dict(reason="warmup_horizon_exhausted")
                        break
                if useful < 3 and warmup_failed is None:
                    warmup_failed = dict(reason="warmup_call_cap")
                row = dict(
                    op=op,
                    useful=useful,
                    calls=calls,
                    state=session.state(),
                    warmup_complete=warmup_failed is None,
                    warmup_failure=warmup_failed,
                )
            elif op == "mark":
                # MacResources.sample walks the complete subtree, including
                # live B12 workers and the spawn helper. Do not add them again.
                row = dict(
                    op=op, marker=session.cpu_marker(), state=session.state()
                )
            elif op == "call" and warmup_failed is not None:
                row = dict(
                    op=op,
                    measured=False,
                    censored=True,
                    reason="warmup_not_admitted",
                    warmup_failure=warmup_failed,
                    state=session.state(),
                )
            elif op == "call":
                row = measured_call(session, request)
                row["measured"] = not bool(row.get("native_assertion_failure"))
            else:
                raise ValueError("unknown finite confirmation operation")
            print(json.dumps(row), flush=True)
            if cold:
                break
    finally:
        session.close()
        close_spawn_helper()


if __name__ == "__main__":
    serve(cold=sys.argv[1:] == ["cold"])
