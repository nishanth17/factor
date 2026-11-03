"""Counterbalanced root and multiplier-setup calibration groups."""

import argparse
import hashlib
import json
import signal
import time

from . import execute_v4, plan_v4


def sample_root(row, fixture, base, multiplier, assignments):
    """Use the production family caller on one shared base/A population."""
    started, cpu = time.monotonic(), time.process_time()
    job, budget, _ = execute_v4.instantiate(
        row,
        fixture,
        execute_v4.runtime_modules(row),
        seconds=row["seconds"],
        work=row["work"],
    )
    job.base, job.multiplier, job.assignments = base, multiplier, assignments
    endpoints = []
    for index in range(8 * 64):
        step = job._next_step()
        if step is None:
            raise RuntimeError("root family exhausted before declared prefix")
        polynomial, roots = step
        if index % 64 == 63:
            endpoints.append((polynomial.a, polynomial.b, len(roots)))
    elapsed = time.monotonic() - started
    consumed = time.process_time() - cpu
    for a, b, _ in endpoints:
        assert (b * b - base.n_prime) % a == 0
    return dict(
        wall=elapsed,
        cpu=consumed,
        work=budget.used,
        fingerprint=hashlib.sha256(str(endpoints).encode()).hexdigest(),
        endpoint_count=len(endpoints),
    )


def sample_setup(row, fixture):
    started, cpu = time.monotonic(), time.process_time()
    job, budget, _ = execute_v4.instantiate(
        row,
        fixture,
        execute_v4.runtime_modules(row),
        seconds=row["seconds"],
        work=row["work"],
    )
    job._setup()
    elapsed = time.monotonic() - started
    consumed = time.process_time() - cpu
    if job.divisor is not None:
        assert 1 < job.divisor < fixture["n"]
        assert fixture["n"] % job.divisor == 0
    scores = job.stats.get("multiplier_scores", ())
    return dict(
        wall=elapsed,
        cpu=consumed,
        work=budget.used,
        k=job.multiplier,
        base_count=0 if job.base is None else len(job.base.entries),
        score_sha256=hashlib.sha256(str(scores).encode()).hexdigest(),
    )


def group(phase, digits):
    plan = json.loads(plan_v4.PLAN.read_text())
    plan_v4.verify_plan(plan)
    rows = [
        r
        for r in plan["calibration"]["calls"]
        if r["phase"] == phase and r["digits"] == digits
    ]
    if len(rows) != (3 if phase == "roots" else 2):
        raise ValueError("group cell census changed")
    corpus = json.loads(plan_v4.CORPUS.read_text())
    fixture = execute_v4.identity(rows[0], corpus)
    context = {}
    slot_work = plan["calibration"]["limits"]["warm_slot_work"]
    for row in rows:
        context[row["id"]] = dict(
            row=row,
            before=execute_v4.warm(
                row,
                fixture,
                execute_v4.runtime_modules(row),
                pre=True,
                work=slot_work,
            ),
            samples=[],
        )
    shared = dict(work=0, wall=0.0, cpu=0.0)
    if phase == "roots":
        first = rows[0]
        common, budget, _ = execute_v4.instantiate(
            first,
            fixture,
            execute_v4.runtime_modules(first),
            seconds=first["seconds"],
            work=first["work"],
        )
        started, cpu = time.monotonic(), time.process_time()
        common._setup()
        if common.divisor is not None:
            raise RuntimeError("shared root setup found exceptional divisor")
        assignments = tuple(common.assignments[index] for index in range(8))
        base, multiplier = common.base, common.multiplier
        shared = dict(
            work=budget.used,
            wall=time.monotonic() - started,
            cpu=time.process_time() - cpu,
        )
    for round_index in range(9):
        if phase == "roots":
            offset = round_index % 3
            order = rows[offset:] + rows[:offset]
        else:
            order = (
                rows if (round_index + (digits == 70)) % 2 == 0 else rows[::-1]
            )
        round_samples = []
        for row in order:
            entry = context[row["id"]]
            if phase == "roots":
                sample = sample_root(
                    row, fixture, base, multiplier, assignments
                )
            else:
                sample = sample_setup(row, fixture)
            entry["samples"].append(sample)
            round_samples.append(sample)
            if sum(x["wall"] for x in entry["samples"]) > row["seconds"]:
                raise RuntimeError("diagnostic wall cap exceeded")
            if sum(x["cpu"] for x in entry["samples"]) > row["seconds"]:
                raise RuntimeError("diagnostic CPU cap exceeded")
            if sum(x["work"] for x in entry["samples"]) > row["work"]:
                raise RuntimeError("diagnostic work cap exceeded")
        if phase == "roots":
            if len({s["fingerprint"] for s in round_samples}) != 1:
                raise RuntimeError("root policies changed physical trials")
    members = []
    for row in rows:
        entry = context[row["id"]]
        after = execute_v4.warm(
            row,
            fixture,
            execute_v4.runtime_modules(row),
            pre=False,
            work=slot_work - entry["before"]["work"],
        )
        measured_work = sum(s["work"] for s in entry["samples"])
        members.append(
            dict(
                id=row["id"],
                policy=(
                    row["configuration"][
                        "root_update_policy"
                        if phase == "roots"
                        else "multiplier_policy"
                    ]
                ),
                before=entry["before"],
                after=after,
                samples=entry["samples"],
                measured_work=measured_work,
                measured_wall=sum(s["wall"] for s in entry["samples"]),
                measured_cpu=sum(s["cpu"] for s in entry["samples"]),
            )
        )
    total_work = shared["work"] + sum(
        m["measured_work"] + m["before"]["work"] + m["after"]["work"]
        for m in members
    )
    return dict(
        phase=phase,
        digits=digits,
        fixture=fixture["id"],
        shared_setup=shared,
        members=members,
        work=total_work,
        reason="diagnostic_complete",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("roots", "ranking"))
    parser.add_argument("digits", type=int, choices=(62, 70))
    args = parser.parse_args()
    count = 3 if args.phase == "roots" else 2
    signal.signal(
        signal.SIGALRM,
        lambda *_: (_ for _ in ()).throw(
            TimeoutError("diagnostic group hard wall allowance")
        ),
    )
    signal.alarm(count * 75 + 5)
    output = group(args.phase, args.digits)
    signal.alarm(0)
    encoded = json.dumps(output, sort_keys=True)
    if len(encoded.encode()) > 64 * 1024:
        raise RuntimeError("recorder_exhaustion: group output")
    print(encoded, flush=True)


if __name__ == "__main__":
    main()
